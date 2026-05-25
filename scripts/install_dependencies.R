#### Install R package dependencies ----
## Run this once after cloning the repo, before any other script.
## Re-runs are cheap: already-installed packages are skipped.

## CRAN packages used anywhere in the codebase (union across all scripts).
## Grouped by purpose; the install order does not matter.
required_packages <- c(
  ## Data wrangling
  "dplyr", "tidyr", "lubridate", "forecast",
  ## Plotting
  "ggplot2", "maps", "mapdata",
  ## Bayesian models (brms also needs a working C++ toolchain via rstan/cmdstanr)
  "brms", "bayestestR",
  ## Geospatial
  "terra", "sf", "rnaturalearth",
  ## Genesys accession API
  "genesysr",
  ## Table rendering (webshot2 also needs Chrome/Chromium installed)
  "kableExtra", "webshot2", "png", "base64enc"
)

missing <- required_packages[!vapply(required_packages, requireNamespace,
                                      logical(1), quietly = TRUE)]

if (length(missing) == 0) {
  message("All required packages are already installed.")
} else {
  message(sprintf("%d of %d required packages are missing:",
                  length(missing), length(required_packages)))
  message(paste0("  - ", missing, collapse = "\n"))
  cat("\nInstall them now from CRAN? (yes/no): ")
  answer <- readline()
  if (tolower(trimws(answer)) %in% c("yes", "y")) {
    install.packages(missing)
    still_missing <- missing[!vapply(missing, requireNamespace,
                                     logical(1), quietly = TRUE)]
    if (length(still_missing) > 0) {
      warning("The following packages failed to install: ",
              paste(still_missing, collapse = ", "),
              ". See the messages above for details.")
    } else {
      message("All packages installed successfully.")
    }
  } else {
    message("Skipping installation. The scripts will fail until these packages ",
            "are installed manually with install.packages().")
  }
}

## Notes on extra system-level prerequisites:
##  * brms needs a working C++ toolchain. On Windows, install Rtools matched to
##    your R version; on macOS, install Xcode command-line tools.
##  * webshot2 needs a Chrome/Chromium browser available on the system path
##    (used by table.R to render HTML tables to PNG).
