# Pest evolution amplifies projected crop losses under climate change

R code and data retrieval accompanying von Schmalensee et al. 2026 "Pest evolution amplifies projected crop losses under climate change"

## Data availability

Raw data, fitted Bayesian models, and the Genesys accession cache are archived on Zenodo at **https://sandbox.zenodo.org/records/502301**. The GitHub repository contains only the code; running `scripts/fetch_data.R` once after cloning will download everything into the appropriate local folders.

## Folder structure

- `data/` — raw experimental and climate data loaded by the scripts.
- `scripts/` — R scripts (model fitting, simulation, figures, tables).
- `models/` — fitted `brms` model objects.
- `output/` — intermediate outputs (parameter estimates, simulation results, the Genesys accession cache) produced by the scripts and retrieved from Zenodo on first run.
- `figures/raw/` — PDF figures produced by the scripts. 

Note: the final figures have been assembled from the raw PDFs and adjusted (aesthetically) in inkscape.

## Running

Scripts are run from the `scripts/` directory and reference paths as `../data/`, `../output/`, etc. `packages_functions_data.R` is sourced by every script and loads shared functions and data.

**First-time setup** — from the `scripts/` directory:

1. `source("install_dependencies.R")` — installs every CRAN package any script will need (installed packages are skipped).
2. `source("fetch_data.R")` — downloads the raw data, fitted models, parameter estimates, simulation outputs, and Genesys cache from Zenodo into `../data/`, `../models/`, and `../output/` (≈3.8 GB total). `fetch_data.R` also works from the repo root; already-present files are skipped on re-runs.

With everything in place, the figure and table scripts run without needing to refit any Bayesian model or rerun any simulation.

With the bundled fitted models, the typical end-to-end workflow is:

1. `temperature_simulation.R` — simulates the 78-year hourly temperature series for the 20 California sites.
2. `results_20_sites.R` and `results_20_sites_daily averages.R` — main and SI simulations of fitness and crop damage across sites.
3. Figure and table scripts (`intro_fig.R`, `empirical_results.R`, `global_maps.R`, `economic_impact.R`, `global_sums.R`, `parameter_table.R`, `table.R`, etc.).

To refit the Bayesian models from scratch, run `development.R`, `growth.R`, `LRS.R` and their `*_global.R` counterparts before the simulation steps; this regenerates everything in `models/` and the parameter-estimate text files in `output/`.

## Dependencies

R (≥ 4.3) plus the following packages:

- **Bayesian models**: `brms`, `bayestestR`
- **Data wrangling**: `dplyr`, `tidyr`, `lubridate`, `forecast`
- **Plotting**: `ggplot2`
- **Geospatial**: `terra`, `sf`, `rnaturalearth`, `maps`, `mapdata`
- **Table rendering**: `kableExtra`, `webshot2`, `png`, `base64enc`
- **Genesys accession API**: `genesysr`

`scripts/install_dependencies.R` installs all of the above in one step. Individual scripts then `library()` only the packages they actually use. `brms` requires a working C++ toolchain (via `rstan` / `cmdstanr`), and `webshot2` requires Chrome/Chromium for the HTML-to-image step used by `table.R`.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The accompanying data and fitted models on Zenodo are released under CC-BY 4.0.
