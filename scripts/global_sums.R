#### Setup ----
## Aggregate global / cowpea-country crop damage under three scenarios:
##   present   = ancestral genotype, present climate
##   warming   = ancestral genotype, future climate
##   evolution = winning genotype (max lambda_geo) per pixel, future climate
## Non-viable pixels (lambda_geo < 1) contribute 0 damage in their scenario.

library(terra)

pre <- rast("../data/global_pres_mean.tif")
fut <- rast("../data/global_fut_mean.tif")

#### Build damage and viability layers ----

dam_anc_pre <- pre[["damage_anc"]]
via_anc_pre <- pre[["lambda_geo_anc"]] >= 1
dam_anc_fut <- fut[["damage_anc"]]
via_anc_fut <- fut[["lambda_geo_anc"]] >= 1

l_anc  <- fut[["lambda_geo_anc"]]
l_hot  <- fut[["lambda_geo_hot"]]
l_cold <- fut[["lambda_geo_cold"]]
l_max  <- max(l_anc, l_hot, l_cold)
is_anc  <- l_anc == l_max
is_hot  <- (l_hot == l_max) & !is_anc
is_cold <- (l_cold == l_max) & !is_anc & !is_hot
via_best_fut <- l_max >= 1
dam_best_fut <- (fut[["damage_anc"]]  * is_anc) +
                (fut[["damage_hot"]]  * is_hot) +
                (fut[["damage_cold"]] * is_cold)

## Mask non-viable pixels to 0 (extinct -> no damage)
mask_zero  <- function(r, v) { x <- r; x[!v] <- 0; x }
dam_pres_v <- mask_zero(dam_anc_pre,  via_anc_pre)
dam_warm_v <- mask_zero(dam_anc_fut,  via_anc_fut)
dam_best_v <- mask_zero(dam_best_fut, via_best_fut)

A <- cellSize(dam_anc_pre, unit = "km") # km^2 per pixel for area-weighting

#### Global means ----
## Use means instead of sums so values are directly comparable to country-level
## means. % changes are invariant under mean vs sum (constant N cancels)

global_summary <- function(P, W, E, scope, aggregation) {
  data.frame(
    scope            = scope,
    aggregation      = aggregation,
    present          = P,
    warming          = W,
    warming_evol     = E,
    pct_warming      = 100 * (W - P) / P,
    pct_warming_evol = 100 * (E - P) / P,
    fold_ratio       = (E - P) / (W - P)
  )
}

## Common pixel mask: only count pixels with defined damage in all three scenarios.
## Ensures a fixed denominator across scenarios.
valid <- !is.na(values(dam_anc_pre)) & !is.na(values(dam_anc_fut)) & !is.na(values(dam_best_fut))
n_valid    <- sum(valid)
A_vals     <- values(A)
total_area <- sum(A_vals[valid])
dp <- values(dam_pres_v); dp[!valid] <- NA
dw <- values(dam_warm_v); dw[!valid] <- NA
de <- values(dam_best_v); de[!valid] <- NA

P_u <- sum(dp, na.rm = TRUE) / n_valid
W_u <- sum(dw, na.rm = TRUE) / n_valid
E_u <- sum(de, na.rm = TRUE) / n_valid
P_a <- sum(dp * A_vals, na.rm = TRUE) / total_area
W_a <- sum(dw * A_vals, na.rm = TRUE) / total_area
E_a <- sum(de * A_vals, na.rm = TRUE) / total_area

global_tab <- rbind(
  global_summary(P_u, W_u, E_u, "Global", "Unweighted grid cell mean"),
  global_summary(P_a, W_a, E_a, "Global", "Area-weighted grid cell mean")
)

#### Cowpea producers (accession-weighted) ----
## Mean damage over the unique accession locations used in Fig 6C:
## producer countries + India retained at n_accessions >= 2 (i.e., the same
## countries that appear in the country_damage CSV with a non-NA country name).
## producer countries + India retained at n_accessions >= 2 (i.e., the same
## countries that appear in the country_damage CSV with a non-NA country name).
## Each unique sampling location contributes once, as in the Fig 6C country means.

cd  <- read.csv("../output/cowpea_country_damage.csv")
acc <- read.csv("../output/cowpea_accession_damage.csv")

fig6c_iso3 <- cd$iso3[!is.na(cd$country)] # producers + India
acc_used <- acc[acc$origcty %in% fig6c_iso3 &
                !is.na(acc$damage_present) &
                !is.na(acc$damage_fut_no_evol) &
                !is.na(acc$damage_fut_with_evol), ]

P_acc <- mean(acc_used$damage_present)
W_acc <- mean(acc_used$damage_fut_no_evol)
E_acc <- mean(acc_used$damage_fut_with_evol)

cowpea_tab <- global_summary(P_acc, W_acc, E_acc,
                             "Cowpea producers", "Accession-weighted mean")

summary_tab <- rbind(global_tab, cowpea_tab)

#### Decomposition: extinction / colonization / magnitude ----
## Decomposes (future_total - present_total) into three additive components.
##   stable      : pixels viable in both periods -- magnitude shift only
##   extinction  : pixels viable in present, non-viable in future (loss = -present)
##   colonization: pixels non-viable in present, viable in future  (gain = +future)
## Sum of components reproduces (W - P) or (E - P).

decompose <- function(dp, vp, df, vf, valid_mask, divisor) {
  ## Restrict viability to common mask; NAs become FALSE so they fall out of all buckets.
  vp <- vp & valid_mask; vp[is.na(vp)] <- FALSE
  vf <- vf & valid_mask; vf[is.na(vf)] <- FALSE
  list(
    stable       = sum(df[vp & vf]   - dp[vp & vf], na.rm = TRUE) / divisor,
    extinction   = sum(-dp[vp & !vf],                na.rm = TRUE) / divisor,
    colonization = sum(df[!vp & vf],                 na.rm = TRUE) / divisor,
    total        = (sum(df, na.rm = TRUE) - sum(dp, na.rm = TRUE)) / divisor
  )
}

vp_anc_pre  <- values(via_anc_pre)
vp_anc_fut  <- values(via_anc_fut)
vp_best_fut <- values(via_best_fut)

dec_warm_u <- decompose(dp, vp_anc_pre, dw, vp_anc_fut,  valid, divisor = n_valid)
dec_evol_u <- decompose(dp, vp_anc_pre, de, vp_best_fut, valid, divisor = n_valid)

## Pixel counts per component (within the common valid mask)
restrict <- function(v) { y <- v & valid; y[is.na(y)] <- FALSE; y }
vp <- restrict(vp_anc_pre)
vw <- restrict(vp_anc_fut)
vb <- restrict(vp_best_fut)

count_components <- function(vp, vf) {
  c(stable       = sum(vp &  vf),
    extinction   = sum(vp & !vf),
    colonization = sum(!vp & vf))
}
n_warm <- count_components(vp, vw)
n_evol <- count_components(vp, vb)

decomp_row <- function(d, n_comp, P, scenario) {
  abs_change <- c(d$stable, d$extinction, d$colonization)
  total      <- d$total
  data.frame(
    scenario        = scenario,
    component       = c("Gradual changes (λ remains ≥ 1)",
                        "Decrease from local extinction (λ goes < 1)",
                        "Increase from local colonizations (λ goes ≥ 1)"),
    n_grid_cells    = n_comp,
    pct_all_pixels  = 100 * n_comp / n_valid,
    absolute_change = abs_change,
    pct_of_present  = 100 * abs_change / P,
    pct_of_total    = 100 * abs_change / total
  )
}

decomp_tab <- rbind(
  decomp_row(dec_warm_u, n_warm, P_u, "Warming only"),
  decomp_row(dec_evol_u, n_evol, P_u, "Warming + evolution")
)

#### Flip effect (warming-driven crop damage decline reversed by evolution) ----
## Counts over the full valid mask (all modelled pixels, n_valid).
## "warm_decreases" = warming alone (ancestor genotype, future vs. present climate)
## yields lower damage than the present baseline. Includes grid cells where warming
## drives the ancestor extinct, since their damage drops to zero (dw = 0 < dp).
## "evol_flips" = within warm_decreases, evolution lifts damage above present.

warm_decreases <- valid & !is.na(dp) & !is.na(dw) & (dw < dp)
evol_flips     <- warm_decreases & !is.na(de) & (de > dp)
flip_tab <- data.frame(
  subset           = c(
    "Grid cells where warming alone decreases damage (ancestor → ancestor)",
    "  ... of which evolution lifts damage above present (flip)"),
  n_grid_cells     = c(sum(warm_decreases), sum(evol_flips)),
  pct_common_mask  = 100 * c(sum(warm_decreases), sum(evol_flips)) / n_valid
)
flip_rate <- 100 * sum(evol_flips) / sum(warm_decreases)

## Flipped grid cells located in Nigeria + Niger (top tier 1 + co-tier 2 producers
## sitting in the Sahel flip zone)
suppressPackageStartupMessages({ library(sf); library(rnaturalearth) })
world      <- ne_countries(scale = "medium", returnclass = "sf")
top_prod   <- vect(world[world$iso_a3 %in% c("NGA", "NER"), ])
top_mask   <- !is.na(values(rasterize(top_prod, dam_anc_pre)))
n_top_flip <- sum(evol_flips & top_mask)
n_top_land <- sum(valid      & top_mask)
flip_top <- list(
  countries       = "Nigeria + Niger",
  n_flip          = n_top_flip,
  pct_global_flip = 100 * n_top_flip / sum(evol_flips),
  pct_global_land = 100 * n_top_land / n_valid,
  enrichment      = (n_top_flip / sum(evol_flips)) / (n_top_land / n_valid)
)
write.csv(as.data.frame(flip_top),
          "../output/global_sums_flip_top_producers.csv", row.names = FALSE)

#### Write summary table ----

dir.create("../output", showWarnings = FALSE)
out_path <- "../output/global_sums_table.txt"

old_width <- options(width = 200); on.exit(options(old_width), add = TRUE)
sink(out_path)
cat("Global crop damage under three scenarios\n")
cat("Generated by scripts/global_sums.R\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("== CROP DAMAGE RATE AND % CHANGE ==\n")
print(summary_tab, row.names = FALSE, digits = 3)

cat("\n== ATTRIBUTION OF CHANGE IN GLOBAL MEAN DAMAGE RATE ==\n")
print(decomp_tab, row.names = FALSE, digits = 3)

cat("\n== SCENARIOS WHERE EVOLUTION FLIPS A WARMING-DRIVEN DAMAGE DECLINE TO AN INCREASE ==\n")
print(flip_tab, row.names = FALSE, digits = 3)
cat(sprintf("\nFlip rate (of grid cells where warming alone decreases damage, %% with evolution-driven damage above present): %.1f%%\n", flip_rate))

cat("\nNotes:\n")
cat("- Present              = ancestral genotype, current climate\n")
cat("- Warming              = ancestral genotype, future climate\n")
cat("- Warming + evolution  = winning genotype (max lambda_geo) per grid cell, future climate\n")
cat("- Grid cells with lambda_geo < 1 contribute 0 damage in their scenario\n")
cat("- Decomposition splits (future - present) into gradual / local extinction / local colonization\n")
cat("- Flip effect: grid cells where warming alone decreases damage but evolution lifts it above present\n")
sink()

## CSV copies
write.csv(summary_tab, "../output/global_sums_aggregate.csv",     row.names = FALSE)
write.csv(decomp_tab,  "../output/global_sums_decomposition.csv", row.names = FALSE)
write.csv(flip_tab,    "../output/global_sums_flip.csv",          row.names = FALSE)

#### Console echo ----

cat("\n== AGGREGATE ==\n")
print(summary_tab, row.names = FALSE, digits = 3)
cat("\n== DECOMPOSITION ==\n")
print(decomp_tab, row.names = FALSE, digits = 3)
cat("\n== FLIP ==\n")
print(flip_tab, row.names = FALSE, digits = 3)
cat(sprintf("Flip rate = %.1f%%\n", flip_rate))
cat("\nWritten to:", out_path, "\n")
cat("CSVs:    ../output/global_sums_aggregate.csv\n")
cat("         ../output/global_sums_decomposition.csv\n")
cat("         ../output/global_sums_flip.csv\n")
