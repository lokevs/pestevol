#### Setup ----

library(brms)
library(bayestestR)

#### Load models ----

devrate <- readRDS("../models/dev_rate_model.RDS")
devrate_global <- readRDS("../models/dev_rate_model_global.RDS")
growthrate <- readRDS("../models/growth_rate_model.RDS")
growthrate_global <- readRDS("../models/growth_rate_model_global.RDS")
offspring <- readRDS("../models/offspring_model.RDS")
offspring_global <- readRDS("../models/offspring_model_global.RDS")

#### Helper functions ----

## Map abbreviated origin/selection names to full names
expand_names <- function(x) {
  map <- c(
    "anc" = "Ancestral",
    "bra" = "Brazil",
    "ca"  = "California",
    "yem" = "Yemen",
    "hot" = "Hot",
    "cold" = "Cold"
  )
  ifelse(x %in% names(map), map[x], x)
}

## Convert scientific notation in prior strings to fixed notation
clean_prior_string <- function(txt_vec) {
  sapply(txt_vec, function(txt) {
    matches <- gregexpr("[0-9]+(\\.[0-9]+)?[eE][+-]?[0-9]+", txt)
    reg_matches <- regmatches(txt, matches)
    if (length(reg_matches[[1]]) == 0) return(txt)

    replacements <- sapply(reg_matches[[1]], function(x) {
      format(as.numeric(x), scientific = FALSE, trim = TRUE)
    })

    res <- txt
    for(i in seq_along(reg_matches[[1]])) {
      res <- sub(reg_matches[[1]][i], replacements[i], res, fixed = TRUE)
    }
    return(res)
  })
}

#### Process model posteriors and priors into a summary table ----

process_model <- function(model_obj, model_pretty_name, model_type) {

  ## Extract priors
  priors_df <- prior_summary(model_obj)
  priors_df$coef[priors_df$coef == ""] <- NA
  priors_df$nlpar[priors_df$nlpar == ""] <- NA

  ## Extract posterior samples (only fixed effects and variance components)
  post <- as.data.frame(model_obj)
  cols_to_keep <- grep("^(b_|sigma|sd_)", names(post), value = TRUE)
  post_subset <- post[, cols_to_keep]

  ## Summarise posteriors with MAP and 90% HDI
  stats <- describe_posterior(
    post_subset,
    centrality = "MAP",
    ci = 0.90,
    ci_method = "HDI"
  )

  stats <- as.data.frame(stats)
  n_rows <- nrow(stats)

  ## Pre-allocate result vectors
  res_model  <- rep(model_pretty_name, n_rows)
  res_type   <- rep(model_type, n_rows)
  res_origin <- character(n_rows)
  res_sel    <- character(n_rows)
  res_param  <- character(n_rows)
  res_prior  <- character(n_rows)

  for (i in 1:n_rows) {
    p <- stats$Parameter[i]

    ## Extract origin metadata
    if (grepl("origin[a-z]+", p)) {
      m <- regexpr("origin[a-z]+", p)
      raw_val <- sub("origin", "", regmatches(p, m))
      res_origin[i] <- expand_names(raw_val)
    } else {
      res_origin[i] <- NA
    }

    ## Extract selection regime metadata
    if (grepl("selection\\.regime[a-z]+", p)) {
      m <- regexpr("selection\\.regime[a-z]+", p)
      raw_val <- sub("selection\\.regime", "", regmatches(p, m))
      res_sel[i] <- expand_names(raw_val)
    } else {
      res_sel[i] <- NA
    }

    ## Clean parameter name
    clean_name <- sub("^b_", "", p)
    clean_name <- sub("replicate__", "", clean_name)
    clean_name <- sub("_Intercept", "", clean_name)
    clean_name <- gsub("origin[a-z]+", "", clean_name)
    clean_name <- gsub("selection\\.regime[a-z]+", "", clean_name)
    clean_name <- gsub(":", "", clean_name)
    clean_name <- sub("_$", "", clean_name)

    if (grepl("^sd_", p)) clean_name <- paste0("SD (", sub("^sd_", "", clean_name), ")")
    if (grepl("^sigma", p)) clean_name <- "Residual SD"

    res_param[i] <- clean_name

    ## Match parameter to its prior
    prior_str <- ""

    if (grepl("^b_", p)) {
      clean_coef <- sub("^b_", "", p)

      match_exact <- priors_df$class == "b" & priors_df$coef == clean_coef
      match_exact[is.na(match_exact)] <- FALSE

      matched_nlpar <- NA
      possible_nlpars <- unique(na.omit(priors_df$nlpar))
      for (nlp in possible_nlpars) {
        if (grepl(paste0("_", nlp, "_"), p) || grepl(paste0("_", nlp, "$"), p)) {
          matched_nlpar <- nlp
          break
        }
      }

      match_nlpar_general <- FALSE
      if (!is.na(matched_nlpar)) {
        match_nlpar_general <- priors_df$class == "b" &
          priors_df$nlpar == matched_nlpar &
          is.na(priors_df$coef)
        match_nlpar_general[is.na(match_nlpar_general)] <- FALSE
      }

      if (any(match_exact)) {
        prior_str <- paste(unique(priors_df$prior[match_exact]), collapse = "; ")
      } else if (any(match_nlpar_general)) {
        prior_str <- paste(unique(priors_df$prior[match_nlpar_general]), collapse = "; ")
      }

    } else if (grepl("^sigma", p)) {
      prior_str <- paste(priors_df$prior[priors_df$class == "sigma"], collapse = "; ")
    } else if (grepl("^sd_", p)) {
      match_sd <- priors_df$class == "sd" & is.na(priors_df$coef)
      match_sd[is.na(match_sd)] <- FALSE
      prior_str <- paste(unique(priors_df$prior[match_sd]), collapse = "; ")
    }

    res_prior[i] <- ifelse(prior_str == "", "Flat", prior_str)
  }

  fmt <- function(x) formatC(x, format = "fg", digits = 3, flag = "#")
  final_priors <- clean_prior_string(res_prior)

  df <- data.frame(
    "Model name" = res_model,
    "Model type" = res_type,
    "Parameter" = res_param,
    "Origin" = res_origin,
    "Selection regime" = res_sel,
    "Mode" = fmt(stats$MAP),
    "Lower 90% HPDI" = fmt(stats$CI_low),
    "Upper 90% HPDI" = fmt(stats$CI_high),
    "Prior" = final_priors,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  return(df)
}

#### Run all models ----

model_list <- list(
  list(devrate, "Development rate", "Origin-specific"),
  list(devrate_global, "Development rate", "Global"),
  list(growthrate, "Growth rate", "Origin-specific"),
  list(growthrate_global, "Growth rate", "Global"),
  list(offspring, "Offspring count", "Origin-specific"),
  list(offspring_global, "Offspring count", "Global")
)

results_table <- do.call(rbind, lapply(model_list, function(x) {
  process_model(x[[1]], x[[2]], x[[3]])
}))

#### Export ----

write.table(results_table, "../output/model_parameters_table.txt", row.names = FALSE, quote = F, sep = "\t")
