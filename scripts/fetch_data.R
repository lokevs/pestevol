#### Fetch data, models, and cached intermediates from Zenodo ----
## Run this once after cloning the repo, before any other script.
## Already-downloaded files are skipped on subsequent runs.
## Works whether sourced from the repo root or from the scripts/ folder.

## Zenodo deposit for von Schmalensee et al. 2026. The concept DOI
## 10.5281/zenodo.20159357 always resolves to the latest version; ZENODO_RECORD
## is the version record that hosts the files.
ZENODO_HOST   <- "zenodo.org"
ZENODO_RECORD <- "20159358"

base_url <- sprintf("https://%s/record/%s/files", ZENODO_HOST, ZENODO_RECORD)

## Raise R's download.file() timeout from its 60 s default to 2 h, so that
## multi-GB files (notably simulated_time_series.txt) don't abort mid-stream
## on residential connections.
options(timeout = max(7200, getOption("timeout")))

## Resolve repo root from whichever directory the user happens to be in
repo_root <- if (basename(normalizePath(getwd(), winslash = "/")) == "scripts") {
  normalizePath("..", winslash = "/")
} else if (file.exists("scripts/fetch_data.R")) {
  normalizePath(".", winslash = "/")
} else {
  stop("Cannot locate the repo root. Run this script with the repo root or ",
       "the scripts/ folder as your working directory.")
}

## Files are uploaded flat to Zenodo; this manifest places each one in the
## correct local folder.
files <- list(
  data = c(
    "thermal_performance.txt",
    "oviposition_rate_LRS.txt",
    "climate_data_california.txt",
    "cowpea_production.txt",
    "global_pres_mean.tif",
    "global_fut_mean.tif",
    "global_tmean_2018_2024.tif",
    "global_tsd_2018_2024.tif"
  ),
  models = c(
    "dev_rate_model.rds",
    "dev_rate_model_global.rds",
    "growth_rate_model.rds",
    "growth_rate_model_global.rds",
    "offspring_model.rds",
    "offspring_model_global.rds"
  ),
  output = c(
    ## Genesys accession cache (skip live API call)
    "cowpea_accessions_geo.csv",
    ## Per-accession + per-country damage tables (used by figure scripts)
    "cowpea_accession_damage.csv",
    "cowpea_country_damage.csv",
    ## Parameter estimates derived from the fitted brms models
    "parameter_estimates_dev_rate.txt",
    "parameter_estimates_dev_rate_global.txt",
    "parameter_estimates_growth_rate.txt",
    "parameter_estimates_growth_rate_global.txt",
    "parameter_estimates_viab_fec.txt",
    "parameter_estimates_viab_fec_global.txt",
    "model_parameters_table.txt",
    ## Table S6 inputs
    "global_sums_aggregate.csv",
    "global_sums_decomposition.csv",
    "global_sums_flip.csv",
    "global_sums_flip_top_producers.csv",
    "global_sums_table.txt",
    ## Simulation outputs (large; skipping these means re-running
    ## temperature_simulation.R + results_20_sites*.R locally)
    "simulated_time_series.txt",
    "validation_time_series.txt",
    "simu_dat_long.txt",
    "simu_dat_long_daily_means.txt"
  )
)

## Ensure all target folders exist (figures/raw/ is needed by figure scripts too)
for (d in c("data", "models", "output", "figures/raw")) {
  dir.create(file.path(repo_root, d), showWarnings = FALSE, recursive = TRUE)
}

## Download every file, skipping any that already exist locally
n_total <- sum(lengths(files))
i <- 0
for (folder in names(files)) {
  for (f in files[[folder]]) {
    i <- i + 1
    dest <- file.path(repo_root, folder, f)
    if (file.exists(dest)) {
      message(sprintf("[%d/%d] %s/%s already present, skipping.", i, n_total, folder, f))
      next
    }
    url <- sprintf("%s/%s?download=1", base_url, f)
    message(sprintf("[%d/%d] Downloading %s -> %s/", i, n_total, f, folder))
    # mode = "wb" is required on Windows for .tif and .rds binaries
    download.file(url, dest, mode = "wb")
  }
}

message("All files retrieved.")
