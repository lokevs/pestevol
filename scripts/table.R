#### Load packages, functions, and data ----
source("packages_functions_data.R")

library("dplyr")
library("tidyr")
library("kableExtra")
library("webshot2")
library("grid")
library("png")
library("base64enc")

#### Compute percentage changes relative to ancestral baseline ----

compute_changes <- function(data, response_var) {
  agg_formula <- as.formula(paste(response_var, "~ trait + selection.regime + origin + site + year + warming.rate"))
  agg <- aggregate(agg_formula, data = data, FUN = "mean")
  agg$change <- 0

  for (i in 1:nrow(agg)) {
    mask <- agg$selection.regime == "anc" &
      agg$site == agg$site[i] &
      agg$origin == agg$origin[i]

    anc_baseline <- agg[[response_var]][mask & agg$warming.rate == 0]
    anc_current  <- agg[[response_var]][mask & agg$warming.rate == agg$warming.rate[i]]

    if (agg$selection.regime[i] == "anc") {
      agg$change[i] <- round(((anc_current / anc_baseline) - 1) * 100)
    } else {
      evo_current <- agg[[response_var]][i]
      agg$change[i] <- round(((evo_current / anc_baseline) - 1) * 100) -
        round(((anc_current / anc_baseline) - 1) * 100)
    }
  }

  agg %>%
    group_by(trait, selection.regime, warming.rate) %>%
    summarise(
      min.value = min(change, na.rm = TRUE),
      lower.quartile = quantile(change, 0.25, na.rm = TRUE),
      mean.value = median(change, na.rm = TRUE),
      upper.quartile = quantile(change, 0.75, na.rm = TRUE),
      max.value = max(change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      selection.regime = case_when(
        selection.regime == "cold" ~ "Cold",
        selection.regime == "anc"  ~ "Ancestral",
        TRUE ~ "Hot"
      ),
      trait = case_when(
        trait == "all"         ~ "All",
        trait == "development" ~ "Development",
        trait == "fecundity"   ~ "Fecundity",
        trait == "growth"      ~ "Growth",
        TRUE                   ~ "Viability"
      )
    )
}

#### Inline boxplot image ----

generate_boxplot_img <- function(min_val, lower_q, mean_val, upper_q, max_val, percent_val, regime) {
  upper_limit <- 2.5  # 250%
  lower_limit <- -1   # -100%

  # Clamp data to visual boundaries
  min_val_lim <- max(min_val/100, lower_limit)
  max_val_lim <- min(max_val/100, upper_limit)

  if (regime == "Hot") colour <- colorRampPalette(c("white", "orangered", "black"))(5)
  if (regime == "Cold") colour <- colorRampPalette(c("white", "dodgerblue", "black"))(5)
  if (regime == "Ancestral") colour <- colorRampPalette(c("white", "darkgrey", "black"))(5)

  f <- tempfile(fileext = ".png")
  png(f, width = 500, height = 200, bg = "white", res = 500)
  par(mar = c(0, 0, 0, 0))

  plot(NULL, xlim = c(lower_limit, upper_limit), ylim = c(-0.1, 1.15), axes = FALSE, xlab = "", ylab = "", type = "n")
  abline(v = 0, col = "grey80", lty = 2)

  # Range rectangle
  rect(min_val_lim, 0.47, max_val_lim, 0.53, col = colour[2], border = colour[2])

  # Overshoot arrows
  if (max_val/100 >= upper_limit) {
    arrows(upper_limit - 0.1, 0.5, upper_limit, 0.5, col = colour[2], lwd = 2, length = 0.1, angle = 20)
  }
  if (min_val/100 <= lower_limit) {
    arrows(lower_limit + 0.1, 0.5, lower_limit, 0.5, col = colour[2], lwd = 2, length = 0.1, angle = 20)
  }

  # IQR rectangle and median line
  lq_lim <- max(lower_q/100, lower_limit)
  uq_lim <- min(upper_q/100, upper_limit)

  rect(lq_lim, 0.3, uq_lim, 0.7, col = colour[2], border = colour[2])
  segments(mean_val/100, 0.3, mean_val/100, 0.7, col = colour[4], lwd = 2)

  # Median percentage label
  text(mean(c(lower_limit, upper_limit) + 0.24), 1.01, font = 2, labels = paste0(ifelse(percent_val > 0, "+", ""), round(percent_val), "%"), pos = 2, cex = 0.75)

  # Range label
  text(upper_limit + 0.23, 0.015, labels = paste0("[", ifelse(round(min_val) > 0, "+", ""), round(min_val), "% to ", ifelse(round(max_val ) > 0, "+", ""), round(max_val), "%]"), pos = 2, cex = 0.65)

  dev.off()

  encoded <- base64enc::dataURI(file = f, mime = "image/png")
  paste0('<img src="', encoded, '" height="40" style="width: auto; max-width: none;"/>')
}

#### Inline summary trend plot ----

generate_summary_plot <- function(plot_specific_data_for_trait, current_selection_regime_for_row, global_ancestral_ref_data) {
  f <- tempfile(fileext = ".png")
  tryCatch({
    png(f, width = 350, height = 150, bg = "white", res = 500)

    current_trait_name <- unique(plot_specific_data_for_trait$trait)[1]
    colors <- c("Ancestral" = "darkgrey", "Hot" = "orangered", "Cold" = "dodgerblue")

    additional_lines_to_draw <- list()

    if (current_trait_name == "All") {
      if (current_selection_regime_for_row == "Hot") {
        data_to_sum <- plot_specific_data_for_trait[plot_specific_data_for_trait$selection.regime == "Hot", ]
        data_to_sum <- data_to_sum[order(data_to_sum$warming.rate), ]
        if (nrow(data_to_sum) > 0 && nrow(global_ancestral_ref_data) > 0 &&
            nrow(data_to_sum) == nrow(global_ancestral_ref_data)) {
          summed_values <- global_ancestral_ref_data$mean.value + data_to_sum$mean.value
          summed_df <- data.frame(warming.rate = global_ancestral_ref_data$warming.rate, mean.value = summed_values)
          additional_lines_to_draw[[length(additional_lines_to_draw) + 1]] <-
            list(data = summed_df, color = alpha(colors["Hot"], 0.7), lty = 3, pch = 16)
        }
      } else if (current_selection_regime_for_row == "Cold") {
        data_to_sum <- plot_specific_data_for_trait[plot_specific_data_for_trait$selection.regime == "Cold", ]
        data_to_sum <- data_to_sum[order(data_to_sum$warming.rate), ]
        if (nrow(data_to_sum) > 0 && nrow(global_ancestral_ref_data) > 0 &&
            nrow(data_to_sum) == nrow(global_ancestral_ref_data)) {
          summed_values <- global_ancestral_ref_data$mean.value + data_to_sum$mean.value
          summed_df <- data.frame(warming.rate = global_ancestral_ref_data$warming.rate, mean.value = summed_values)
          additional_lines_to_draw[[length(additional_lines_to_draw) + 1]] <-
            list(data = summed_df, color = alpha(colors["Cold"], 0.7), lty = 3, pch = 16)
        }
      }
    } else {
      if (current_selection_regime_for_row != "Ancestral") {
        data_to_sum <- plot_specific_data_for_trait[plot_specific_data_for_trait$selection.regime == current_selection_regime_for_row, ]
        data_to_sum <- data_to_sum[order(data_to_sum$warming.rate), ]
        if (nrow(data_to_sum) > 0 && nrow(global_ancestral_ref_data) > 0 &&
            nrow(data_to_sum) == nrow(global_ancestral_ref_data)) {
          summed_values <- global_ancestral_ref_data$mean.value + data_to_sum$mean.value
          summed_df <- data.frame(warming.rate = global_ancestral_ref_data$warming.rate, mean.value = summed_values)
          additional_lines_to_draw[[length(additional_lines_to_draw) + 1]] <-
            list(data = summed_df, color = alpha(colors[current_selection_regime_for_row], 0.7), lty = 3, pch = 16)
        }
      }
    }

    par(mar = c(0, 0, 0, 0))
    plot(NULL, xlim = c(-0.003, 0.063), ylim = c(-12, 250), axes = FALSE, xlab = "", ylab = "")

    if (nrow(global_ancestral_ref_data) > 0) {
      lines(global_ancestral_ref_data$warming.rate, global_ancestral_ref_data$mean.value,
            col = colors["Ancestral"], lwd = 1.7, lty = 1)
      points(global_ancestral_ref_data$warming.rate, global_ancestral_ref_data$mean.value,
             col = colors["Ancestral"], pch = 16, cex = 0.5)
    }

    if (length(additional_lines_to_draw) > 0) {
      for (line_info in additional_lines_to_draw) {
        if (nrow(line_info$data) > 0) {
          lines(line_info$data$warming.rate, line_info$data$mean.value,
                col = line_info$color, lwd = 1.7, lty = line_info$lty)
          points(line_info$data$warming.rate, line_info$data$mean.value,
                 col = line_info$color, pch = line_info$pch, cex = 0.5)
        }
      }
    }
    dev.off()
  }, error = function(e) {
    if (exists("f") && !is.null(f) && !is.null(dev.list()) && (dev.cur() != 1)) {
      suppressWarnings(dev.off())
    }
    stop(e)
  }, finally = {})

  encoded <- base64enc::dataURI(file = f, mime = "image/png")
  unlink(f)
  paste0('<img src="', encoded, '" height="40" style="width: auto; max-width: none;"/>')
}

#### Build boxplot table (full or "All"-only mini) ----

# trait_filter = NULL  -> full 11-row table (1 Anc + 5 Hot + 5 Cold), with rotated row group labels
# trait_filter = "All" -> 3-row mini table (Anc/Hot/Cold for "All" trait only), no rotated labels
build_boxplot_table <- function(agg_means, output_file, trait_filter = NULL) {
  full_table <- is.null(trait_filter)
  if (!full_table) {
    agg_means <- filter(agg_means, trait %in% trait_filter)
  }

  ## Global ancestral reference for summary plots
  global_anc <- agg_means %>%
    filter(trait == "All", selection.regime == "Ancestral") %>%
    arrange(warming.rate)

  ## Summary plots per trait x regime
  base_data <- agg_means %>%
    mutate(selection.regime = factor(selection.regime, levels = c("Ancestral", "Hot", "Cold")))

  summary_col <- base_data %>%
    group_by(trait, selection.regime) %>%
    group_map(~{
      all_data <- base_data[base_data$trait == .y$trait, ]
      tibble(
        trait = .y$trait,
        selection.regime = .y$selection.regime,
        summary_plot_html = generate_summary_plot(all_data, .y$selection.regime, global_anc)
      )
    }, .keep = TRUE) %>%
    bind_rows()

  ## Pivot boxplots into warming-rate columns
  plot_data <- agg_means %>%
    mutate(warming = case_when(
      warming.rate == 0    ~ "0°C/y",
      warming.rate == 0.02 ~ "0.02°C/y",
      warming.rate == 0.04 ~ "0.04°C/y",
      warming.rate == 0.06 ~ "0.06°C/y"
    )) %>%
    rowwise() %>%
    mutate(cell_plot = generate_boxplot_img(
      min.value, lower.quartile, mean.value, upper.quartile, max.value,
      mean.value, selection.regime
    )) %>%
    ungroup() %>%
    select(selection.regime, trait, warming, cell_plot) %>%
    pivot_wider(names_from = warming, values_from = cell_plot) %>%
    mutate(selection.regime = factor(selection.regime, levels = c("Ancestral", "Hot", "Cold"))) %>%
    arrange(selection.regime, trait) %>%
    left_join(summary_col, by = c("trait", "selection.regime"))

  if (full_table) {
    ## Full 11-row table: rotated row group labels, group separators at rows 1, 6, 11
    final <- plot_data %>%
      group_by(selection.regime) %>%
      mutate(
        Group = case_when(
          selection.regime == "Ancestral" & row_number() == 1 ~
            "<div style='writing-mode: vertical-rl; transform: rotate(180deg); font-size: 16px;'>Baseline</div>",
          selection.regime == "Cold" & row_number() == 1 ~
            "<div id='evo-label' style='position: absolute; writing-mode: vertical-rl; transform: rotate(180deg); font-size: 16px;'>Evolutionary deviations from baseline</div>",
          TRUE ~ ""
        ),
        `Selection regime` = ifelse(row_number() == 1, as.character(selection.regime), "")
      ) %>%
      ungroup() %>%
      select(Group, `Selection regime`, Trait = trait,
             `0°C/y`, `0.02°C/y`, `0.04°C/y`, `0.06°C/y`,
             `Summary Plot` = summary_plot_html)

    kbl_out <- final %>%
      kbl(escape = FALSE, align = "lllccccc", format = "html",
          col.names = c("", "Selection<br>regime", "Trait",
                         "0°C/y", "0.02°C/y", "0.04°C/y", "0.06°C/y", "")) %>%
      add_header_above(c(" " = 3, "Warming rate" = 4, " " = 1)) %>%
      kable_styling(full_width = FALSE, bootstrap_options = "hover", position = "left") %>%
      row_spec(0:nrow(plot_data), extra_css = "border: none;") %>%
      column_spec(4:6, extra_css = "border-right: 1px solid lightgrey;") %>%
      column_spec(c(1, 3, 7, 8), extra_css = "border-right: 1px solid black;")

    kbl_out <- gsub("</table>",
      "<style>
      .table thead th { border-bottom: none !important; }
      .table thead th:nth-child(n+2):nth-child(-n+8) { border-bottom: 1px solid black !important; }
      .table tbody tr td { height: 62px !important; line-height: 40px !important; padding-top: 0 !important; padding-bottom: 0 !important; vertical-align: middle !important; white-space: nowrap !important; overflow: visible !important; }
      .table td:first-child { position: relative !important; }
      #evo-label { position: absolute !important; z-index: 10 !important; left: 10px !important; height: 125px !important; }
      .table tbody tr:nth-child(n+2):nth-child(-n+11) td:nth-child(-n+2) { border-bottom: none !important; }
      .table tbody tr:nth-child(n+2):nth-child(-n+11) td:nth-child(n+3):nth-child(-n+8) { border-bottom: 1px solid lightgrey !important; }
      .table tbody tr:nth-child(n+1):nth-child(-n+1) td:nth-child(n+2):nth-child(-n+8),
      .table tbody tr:nth-child(n+6):nth-child(-n+6) td:nth-child(n+2):nth-child(-n+8),
      .table tbody tr:nth-child(n+11):nth-child(-n+11) td:nth-child(n+2):nth-child(-n+8) { border-bottom: 1px solid black !important; }
      </style></table>",
      kbl_out)
  } else {
    ## Mini table: 3 rows (one per regime, "All" trait), every row gets a black border
    final <- plot_data %>%
      mutate(`Selection regime` = as.character(selection.regime)) %>%
      select(`Selection regime`, Trait = trait,
             `0°C/y`, `0.02°C/y`, `0.04°C/y`, `0.06°C/y`,
             `Summary Plot` = summary_plot_html)

    kbl_out <- final %>%
      kbl(escape = FALSE, align = "llccccc", format = "html",
          col.names = c("Selection<br>regime", "Trait",
                         "0°C/y", "0.02°C/y", "0.04°C/y", "0.06°C/y", "")) %>%
      add_header_above(c(" " = 2, "Warming rate" = 4, " " = 1)) %>%
      kable_styling(full_width = FALSE, bootstrap_options = "hover", position = "left") %>%
      row_spec(0:nrow(plot_data), extra_css = "border: none;") %>%
      column_spec(3:5, extra_css = "border-right: 1px solid lightgrey;") %>%
      column_spec(c(2, 6, 7), extra_css = "border-right: 1px solid black;")

    kbl_out <- gsub("</table>",
      "<style>
      .table thead th { border-bottom: none !important; }
      .table thead th:nth-child(n+1):nth-child(-n+7) { border-bottom: 1px solid black !important; }
      .table tbody tr td { height: 62px !important; line-height: 40px !important; padding-top: 0 !important; padding-bottom: 0 !important; vertical-align: middle !important; white-space: nowrap !important; overflow: visible !important; border-bottom: 1px solid black !important; }
      </style></table>",
      kbl_out)
  }

  save_kable(kbl_out, output_file)
  cat("  Saved:", output_file, "\n")
}

#### Load simulation data ----

cat("Loading hourly simulation data (this may take a moment)...\n")
simu_hourly <- read.delim("../output/simu_dat_long.txt")

cat("Loading daily-average simulation data...\n")
simu_daily <- read.delim("../output/simu_dat_long_daily_means.txt")

datasets <- list(
  med   = subset(simu_hourly, mortality == "med"),
  high  = subset(simu_hourly, mortality == "high"),
  low   = subset(simu_hourly, mortality == "low"),
  daily = subset(simu_daily,  mortality == "med")
)

cat("Computing changes for each (dataset, metric)...\n")
crop <- lapply(datasets, compute_changes, response_var = "offspring.growth.rate")
fit  <- lapply(datasets, compute_changes, response_var = "fitness.rate")

#### Tables S1, S2: full boxplot tables (hourly, intermediate mortality) ----

cat("Generating Table S1 (crop damage, full)...\n")
build_boxplot_table(crop$med, "../figures/raw/tab_S1_crop_damage.pdf")

cat("Generating Table S2 (fitness, full)...\n")
build_boxplot_table(fit$med, "../figures/raw/tab_S2_fitness.pdf")

#### Mini boxplot tables for Tables S3 and S4 ----
# Each mini table has 3 rows (Anc/Hot/Cold for the "All" trait) and is a layout
# component to be assembled in Inkscape.
#
# Layout map:
#   Table S3 (extrinsic mortality):
#     intermediate -> tab_mini_crop_med    + tab_mini_fitness_med
#     high         -> tab_mini_crop_high   + tab_mini_fitness_high
#     low          -> tab_mini_crop_low    + tab_mini_fitness_low
#   Table S4 (hourly vs daily-average temperatures):
#     hourly       -> tab_mini_crop_med    + tab_mini_fitness_med   (re-used from S3)
#     daily        -> tab_mini_crop_daily  + tab_mini_fitness_daily

cat("Generating mini tables for S3 + S4...\n")
for (k in names(datasets)) {
  build_boxplot_table(crop[[k]],
                      sprintf("../figures/raw/tab_mini_crop_%s.pdf",    k),
                      trait_filter = "All")
  build_boxplot_table(fit[[k]],
                      sprintf("../figures/raw/tab_mini_fitness_%s.pdf", k),
                      trait_filter = "All")
}

cat("All tables generated.\n")
