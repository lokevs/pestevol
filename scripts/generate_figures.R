#### Generate all figures ----
# Produces all PDFs in figures/raw/, all parameter-estimate files in output/,
# and all table outputs. Sources each script in dependency order.

# Set working directory to script location
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

cat("=== Generating all figures ===\n\n")

run_script <- function(name) {
  cat(sprintf("--- Running %s ---\n", name))
  tryCatch(source(name),
           error = function(e) cat(sprintf("  ERROR: %s\n", e$message)))
  cat(sprintf("  Done.\n\n"))
}

#### Load shared requirements ----
cat("Loading packages and shared data...\n")
source("packages_functions_data.R")
cat("Done loading shared data.\n\n")

#### Model scripts (use cached models in ../models/; refit if missing) ----
## Each writes parameter_estimates_*.txt to ../output/ and PPC/posterior PDFs to ../figures/raw/.
run_script("development.R")
run_script("development_global.R")
run_script("growth.R")
run_script("growth_global.R")
run_script("LRS.R")
run_script("LRS_global.R")

#### Temperature simulation (generates ../output/simulated_time_series.txt) ----
run_script("temperature_simulation.R")

## Pre-load the large simulated time series so downstream scripts reuse the
## in-memory copy instead of re-reading the ~2 GB file.
if (!exists("simulated_time_series")) {
  cat("Loading simulated time series (this takes a minute)...\n")
  simulated_time_series  <- read.delim("../output/simulated_time_series.txt")
  simulated_temperatures <- simulated_time_series # alias used by some scripts
}

#### Main simulation scripts (write simu_dat_long*.txt to ../output/) ----
run_script("results_20_sites.R")
run_script("results_20_sites_daily averages.R")

#### Figure scripts ----
run_script("intro_fig.R")
run_script("empirical_results.R")
run_script("trait_specific.R")
run_script("hourly_rates_global.R")
run_script("global_maps.R")
run_script("economic_impact.R") # uses ../output/cowpea_accessions_geo.csv cache; queries Genesys if missing
run_script("methods_plot.R")

#### Table scripts ----
run_script("table.R")
run_script("parameter_table.R")
run_script("global_sums.R")

cat("=== All outputs generated. Check figures/raw/ and output/ ===\n")
