#### Load packages ----
## Installation is handled separately by install_dependencies.R; this script
## only attaches the packages it (and downstream scripts that source it) need.
library(dplyr)
library(tidyr)
library(lubridate)
library(forecast)
library(ggplot2)
library(brms)
library(bayestestR)  # point_estimate() for MAP estimates in the model scripts

#### Temperature simulation function ----
simulate_temperature <- function(input_temperature, input_date_time, input_site, years = 1, n = 1, warming_rate = 0, start_day = 1, end_day = 365) {
  if(length(warming_rate) == 1) {
    warming_rate <- c(warming_rate)
  }
  logK <- log(input_temperature + 273.15)

  tempdat <- data.frame(temp = logK, date_time = input_date_time, site = input_site)
  tempdat$hour <- hour(tempdat$date_time)
  tempdat$day_of_year <- yday(tempdat$date_time)
  tempdat$block <- ceiling(tempdat$day_of_year / 5)

  ## Block-based hourly variations and daily fluctuations
  get_hourly_daily_variations <- function(data) {
    hourly_avg_by_block <- data %>%
      group_by(block, hour) %>%
      summarise(logK = mean(temp), .groups = 'drop')

    hourly_avg_by_block <- hourly_avg_by_block %>%
      group_by(block) %>%
      mutate(logK = logK - mean(logK)) %>%
      ungroup()

    daily_avg <- data %>%
      group_by(day_of_year) %>%
      summarise(logK = mean(temp), .groups = 'drop')

    daily_avg <- rbind(daily_avg, data.frame(day_of_year = daily_avg$day_of_year + 365, logK = daily_avg$logK))
    daily_spline <- smooth.spline(daily_avg$day_of_year, daily_avg$logK, spar = 0.7)

    list(hourly_avg_by_block = hourly_avg_by_block, daily_spline = daily_spline)
  }

  ## Build day sequence across years
  daylist <- function(start_day, end_day, years){
    if(end_day <= start_day){
      return(c(start_day:365, rep(1:365, years-1), 1:end_day))
    }
    else if(years == 1){
      return(start_day:end_day)
    }
    else{
      return(c(start_day:365, rep(1:365, years-2), 1:end_day))
    }
  }

  ## Simulate ARIMA residuals
  simulate_residuals <- function(residual, length_out) {
    fit <- arima(residual, order = c(1, 0, 0), seasonal = list(order = c(1, 0, 0), period = 24))
    arima_sim <- simulate(fit, nsim = length_out)
    arima_sim <- arima_sim * (sd(residual) / sd(arima_sim))
    return(arima_sim)
  }

  simulations_list <- vector("list", n)
  unique_sites <- unique(input_site)

  for (i in 1:n) {
    if (i > length(unique_sites)) {
      warning(paste("Requested n =", n, "simulations, but only", length(unique_sites), "unique sites available. Stopping at site", i-1))
      simulations_list <- simulations_list[1:(i-1)]
      break
    }

    sampled_site <- unique_sites[i]
    site_data <- tempdat %>% filter(site == sampled_site)
    variations <- get_hourly_daily_variations(site_data)

    time_grid <- expand.grid(hour = 0:23, day = daylist(start_day, end_day, years))
    time_grid$block <- ceiling(time_grid$day %% 365 / 5)
    time_grid$block[time_grid$block == 0] <- 73 # Handle day 365

    daily_vals <- predict(variations$daily_spline, time_grid$day)$y

    # Look up hourly variation for each block
    hourly_vals <- numeric(nrow(time_grid))
    for (j in 1:nrow(time_grid)) {
      current_block <- time_grid$block[j]
      current_hour <- time_grid$hour[j]

      matching_entry <- variations$hourly_avg_by_block %>%
        filter(block == current_block, hour == current_hour)

      if (nrow(matching_entry) > 0) {
        hourly_vals[j] <- matching_entry$logK
      } else {
        closest_block <- variations$hourly_avg_by_block %>%
          filter(hour == current_hour) %>%
          arrange(abs(block - current_block)) %>%
          slice(1)

        if (nrow(closest_block) > 0) {
          hourly_vals[j] <- closest_block$logK
        } else {
          hourly_vals[j] <- 0
        }
      }
    }

    combined_logK <- daily_vals + hourly_vals

    # Calculate residuals using block-specific hourly averages
    site_data$block <- ceiling(site_data$day_of_year / 5)
    residuals <- numeric(nrow(site_data))
    for (j in 1:nrow(site_data)) {
      current_block <- site_data$block[j]
      current_hour <- site_data$hour[j]

      hourly_val <- variations$hourly_avg_by_block %>%
        filter(block == current_block, hour == current_hour)

      hourly_pred <- if (nrow(hourly_val) > 0) hourly_val$logK else 0
      daily_pred <- predict(variations$daily_spline, site_data$day_of_year[j])$y
      residuals[j] <- site_data$temp[j] - (daily_pred + hourly_pred)
    }

    simulated_residuals <- simulate_residuals(residuals, length(combined_logK))
    simulated_logK <- combined_logK + as.numeric(simulated_residuals)
    simulated_temp <- exp(simulated_logK) - 273.15
    base_simulated_temp <- simulated_temp

    simulation_time <- time_grid
    simulation_time$year <- rep(1:years, each = 365 * 24)

    sim_frames <- vector("list", length(warming_rate))

    for (w in 1:length(warming_rate)) {
      warming_per_hour <- warming_rate[w] / (365 * 24)
      current_simulated_temp <- base_simulated_temp + (0:(length(base_simulated_temp) - 1)) * warming_per_hour

      sim_frames[[w]] <- data.frame(
        time = 1:length(current_simulated_temp),
        temperature = as.numeric(current_simulated_temp),
        site = sampled_site,
        simulation = i,
        year = simulation_time$year,
        day = simulation_time$day + (simulation_time$year - 1) * 365,
        hour = simulation_time$hour,
        yday = simulation_time$day,
        block = simulation_time$block,
        warming.rate = warming_rate[w]
      )
    }

    simulations_list[[i]] <- do.call(rbind, sim_frames)
    print(paste0(i, " of ", n))
  }

  simulations <- do.call(rbind, simulations_list)
  return(simulations)
}

#### Helper functions ----

geomean <- function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

invlogit <- function(x)(1/(1+exp(-x)))

#### Thermal performance curve functions ----

## Lobry–Rosso–Flandrois (LRF) development/growth rate function
LRF <- function(temp, Tmin, Tmax, Topt, Ropt){
  ifelse(temp > Tmax | temp < Tmin, 0, Ropt * ((temp + 273.15) - (Tmax + 273.15)) * ((temp + 273.15) - (Tmin + 273.15)) ^ 2 / (((Topt + 273.15) - (Tmin + 273.15)) * (((Topt + 273.15) - (Tmin + 273.15)) * ((temp + 273.15) - (Topt + 273.15)) - ((Topt + 273.15) - (Tmax + 273.15)) * ((Topt + 273.15) + (Tmin + 273.15) - 2 *(temp + 273.15)))))
}
LRF <- Vectorize(LRF)

## Fecundity function (lifetime reproductive success)
fec_function <- function(temp, Topt, fecmax, breadth, fatness) {
  20 ^ (-abs((temp - Topt) / (breadth/2)) ^ (fatness + 1)) * fecmax
}
fec_function <- Vectorize(fec_function)

## Viability function (logit-scale mortality)
viab_function <- function(temp, a, h, k){
  a * (temp - h)^2 + k
}
viab_function <- Vectorize(viab_function)

## Viability rate (daily)
viability_rate <- function(temp, Tmin, Tmax, Topt, Ropt, a, h, k){
  dev.rate = ifelse(temp > 37, LRF(37, Tmin, Tmax, Topt, Ropt),
                    ifelse(temp < 17, LRF(17, Tmin, Tmax, Topt, Ropt),
                           LRF(temp, Tmin, Tmax, Topt, Ropt)))
  prop.viable = 1 - invlogit(viab_function(temp, a, h, k))
  return(prop.viable^dev.rate)
}

#### Stan functions for brms models ----
stan_funs <- "
  real fec_function(real temp, real Topt, real fecmax, real breadth, real fatness) {
    return 20 ^ (-abs((temp - Topt) / (breadth/2)) ^ (fatness + 1)) * fecmax;
  }

  real LRF(real temp, real Tmin, real Tmax, real Topt, real Ropt) {
    real devrate;
    if(temp > Tmax||temp < Tmin || Topt - Tmin < Tmax - Topt || Tmin >= Topt || Tmax <= Topt) {
      devrate = 0;
    }
    else{
      devrate = ((Ropt * ((Tmax + 273.15) - (temp + 273.15)) * ((temp + 273.15) - (Tmin + 273.15)) ^ 2) / (((Tmin + 273.15) - (Topt + 273.15)) * ( - (temp + 273.15) * (Tmin + 273.15) + 3 * (temp + 273.15) * (Topt + 273.15) - 2 * (Topt + 273.15) ^ 2 + (Tmax + 273.15) * ( - 2 * (temp + 273.15) + (Tmin + 273.15) + (Topt + 273.15)))));
    }
    return devrate;
  }

  real hurdle_gaussian_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             normal_lpdf(y | mu, sigma);
    }
  }
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

#### Read and prepare data ----
dat <- read.delim("../data/thermal_performance.txt")

## Dummy variables for model fitting
dat$random.dummy <- ifelse(dat$replicate == "anc.yem" | dat$replicate == "anc.ca", 0, 1)
dat$anc.dummy  <- as.factor(ifelse(dat$selection.regime == "anc", "anc", "not"))

## Derived traits
dat$dev.rate <- 1/as.numeric(as.POSIXct(dat$first.hatch.date)-as.POSIXct(dat$start.date))
dat$temperature <- as.numeric(dat$temperature)
dat$mean.weight <- dat$weight/dat$adult.offspring
dat$growth.rate <- dat$mean.weight * dat$dev.rate
dat$adult.offspring <- dat$adult.offspring - 2 # remove the parents
dat$temperature.factor <- as.factor(dat$temperature)

## Plotting variables
dat$colour <- NA
dat$point <- NA
dat$linetype <- NA
for(i in 1:nrow(dat)){
  if(dat$origin[i] == "bra") { dat$point[i] <- 21; dat$linetype[i] <- 1 }
  if(dat$origin[i] == "yem") { dat$point[i] <- 22; dat$linetype[i] <- 2 }
  if(dat$origin[i] == "ca")  { dat$point[i] <- 24; dat$linetype[i] <- 3 }
  if(dat$selection.regime[i] == "hot")      dat$colour[i] <- "orangered"
  if(dat$selection.regime[i] == "cold")     dat$colour[i] <- "dodgerblue"
  if(dat$selection.regime[i] == "anc")      dat$colour[i] <- "darkgrey"
}

dat$relative.adult.offspring <- 0
for(i in 1:nrow(dat)){
  dat$relative.adult.offspring[i] <- dat$adult.offspring[i] / mean(dat$adult.offspring[dat$selection.regime == dat$selection.regime[i] & dat$temperature == dat$temperature[i] & dat$origin == dat$origin[i]], na.rm = T)
}

#### Oviposition rate / LRS data ----
oviposition_rate_LRS <- read.delim("../data/oviposition_rate_LRS.txt")
oviposition_rate_LRS$ratio <- oviposition_rate_LRS$eggs.1h / oviposition_rate_LRS$LRS
oviposition_rate_LRS$log.ratio <- log(oviposition_rate_LRS$ratio)

# Batch correction using replicated Yemen switch line
cor_fm1 <- lm(log.ratio ~ I(temperature - 29), subset(oviposition_rate_LRS, selection.regime == "sw90"))
cor_fm2 <- lm(log.ratio ~ I(temperature - 29), subset(oviposition_rate_LRS, selection.regime == "sw"))
intercept_shift <- coef(cor_fm2)[1] - coef(cor_fm1)[1]
oviposition_rate_LRS$log.ratio.corrected <- ifelse(oviposition_rate_LRS$selection.regime == "sw90" | oviposition_rate_LRS$selection.regime == "anc",
                                         oviposition_rate_LRS$log.ratio + intercept_shift,
                                         oviposition_rate_LRS$log.ratio)

oviposition_rate_LRS <- subset(oviposition_rate_LRS, selection.regime != "sw" & selection.regime != "sw90")

#### Oviposition rate functions ----

## Daily rates
oviposition_rate_function_anc <- function(temp, LRS){
  coefs <- coef(lm(log.ratio.corrected ~ temperature, subset(oviposition_rate_LRS, selection.regime == "anc")))
  return(exp(coefs[1] + coefs[2] * temp) * LRS * 24)
}
oviposition_rate_function_selection <- function(temp, LRS){
  coefs <- coef(lm(log.ratio.corrected ~ temperature, subset(oviposition_rate_LRS, selection.regime != "anc")))
  return(exp(coefs[1] + coefs[2] * temp) * LRS* 24)
}

## Hourly rates (pre-baked coefficients for speed)
coefs_anc <- coef(lm(log.ratio.corrected ~ temperature, subset(oviposition_rate_LRS, selection.regime == "anc")))
oviposition_rate_function_anc_hourly <- eval(parse(text = paste0('function(temp, LRS){
  return(exp(', coefs_anc[1], ' + ', coefs_anc[2], ' * temp) * LRS)
}')))

coefs_selection <- coef(lm(log.ratio.corrected ~ temperature, subset(oviposition_rate_LRS, selection.regime != "anc")))
oviposition_rate_function_selection_hourly <- eval(parse(text = paste0('function(temp, LRS){
  return(exp(', coefs_selection[1], ' + ', coefs_selection[2], ' * temp) * LRS)
}')))

#### Hourly rate functions ----

LRF_hourly <- function(temp, Tmin, Tmax, Topt, Ropt){
  ifelse(temp > Tmax | temp < Tmin, 0, (Ropt / 24) * ((temp + 273.15) - (Tmax + 273.15)) * ((temp + 273.15) - (Tmin + 273.15)) ^ 2 / (((Topt + 273.15) - (Tmin + 273.15)) * (((Topt + 273.15) - (Tmin + 273.15)) * ((temp + 273.15) - (Topt + 273.15)) - ((Topt + 273.15) - (Tmax + 273.15)) * ((Topt + 273.15) + (Tmin + 273.15) - 2 *(temp + 273.15)))))
}

viability_rate_hourly <- function(temp, Tmin, Tmax, Topt, Ropt, a, h, k){
  dev.rate = ifelse(temp > 37, LRF_hourly(37, Tmin, Tmax, Topt, Ropt),
                    ifelse(temp < 17, LRF_hourly(17, Tmin, Tmax, Topt, Ropt),
                           LRF_hourly(temp, Tmin, Tmax, Topt, Ropt)))
  prop.viable = 1 - invlogit(viab_function(temp, a, h, k))
  return(prop.viable^dev.rate)
}

#### Fecundity-mass scaling model ----
mass_offspring <- lm(relative.adult.offspring ~ log(mean.weight) + selection.regime, data = subset(dat, temperature == 29))

# # Plot for SI
# pdf("../figures/raw/fig_S19_mass_scaling.pdf", height = 2, width = 1.8, pointsize = 3)
# plot(relative.adult.offspring~mean.weight, data = subset(dat, temperature == 29), pch = 21, ylim = c(0.35, 1.55), xlim = c(0.001, 0.003), bg = colour,ylab = "Proportion of average max fecundity", xlab  = "Mean dry mass (g)")
# regimes <- unique(dat$selection.regime)
# for (i in seq_along(regimes)) {
#   cols <- ifelse(regimes[i] == "anc", "darkgrey", ifelse(regimes[i] == "hot", "orangered", "dodgerblue"))
#   curve(
#     predict(mass_offspring,
#             newdata = data.frame(mean.weight = x,
#                                  selection.regime = regimes[i])),
#     from = 0.001, to = 0.003,
#     add  = TRUE,
#     col  = cols,
#     lwd  = 2
#   )
# }
# dev.off()

#### TPC generator ----
# Builds rate functions from fitted parameter estimates and assigns them to global environment
generate_TPCs <- function(line_origin, regime_fec, regime_viability, regime_dev, regime_growth){

  parameter_subset <- subset(TPC_parameters, selection.regime == regime_fec & trait == "fecundity" & origin == line_origin|
                               selection.regime == regime_viability & trait == "viability" & origin == line_origin |
                               selection.regime == regime_dev & trait == "development_rate" & origin == line_origin |
                               selection.regime == regime_growth & trait == "growth_rate" & origin == line_origin)

  fecmax <<- parameter_subset$MAP[parameter_subset$trait == "fecundity" & parameter_subset$parameter == "fecmax"]

  fecmax_scaling <<- eval(parse(text = paste0('function(x) {',ifelse(regime_fec == "anc", coef(mass_offspring)[1],
                                                    ifelse(regime_fec == "cold", coef(mass_offspring)[1] + coef(mass_offspring)[3],
                                                           coef(mass_offspring)[1] + coef(mass_offspring)[4])), ' + ', coef(mass_offspring)[2], ' * log(x)}')))

  ovi_rate <<-  eval(parse(text = paste0('function(x) {
  LRS <- fec_function(x,
                      Topt = ', parameter_subset$MAP[parameter_subset$trait == "fecundity" & parameter_subset$parameter == "Topt"], ',
                      fecmax = ', parameter_subset$MAP[parameter_subset$trait == "fecundity" & parameter_subset$parameter == "fecmax"], ',
                      breadth = ', parameter_subset$MAP[parameter_subset$trait == "fecundity" & parameter_subset$parameter == "breadth"], ',
                      fatness = ', parameter_subset$MAP[parameter_subset$trait == "fecundity" & parameter_subset$parameter == "fatness"], ')
  if (regime_fec == "anc") {
    return(oviposition_rate_function_anc_hourly(x, LRS))
  } else {
    return(oviposition_rate_function_selection_hourly(x, LRS))
  }
}')))

  dev_rate <<- eval(parse(text = paste0('function(x) {
  return(LRF_hourly(x,
      Tmin = ', parameter_subset$MAP[parameter_subset$trait == "development_rate" & parameter_subset$parameter == "Tmin"], ',
      Tmax = ', parameter_subset$MAP[parameter_subset$trait == "development_rate" & parameter_subset$parameter == "Tmax"], ',
      Topt = ', parameter_subset$MAP[parameter_subset$trait == "development_rate" & parameter_subset$parameter == "Topt"], ',
      Ropt = ', parameter_subset$MAP[parameter_subset$trait == "development_rate" & parameter_subset$parameter == "Ropt"], '))
}')))

  growth_rate <<- eval(parse(text = paste0('function(x) {
  return(LRF_hourly(x,
      Tmin = ', parameter_subset$MAP[parameter_subset$trait == "growth_rate" & parameter_subset$parameter == "Tmin"], ',
      Tmax = ', parameter_subset$MAP[parameter_subset$trait == "growth_rate" & parameter_subset$parameter == "Tmax"], ',
      Topt = ', parameter_subset$MAP[parameter_subset$trait == "growth_rate" & parameter_subset$parameter == "Topt"], ',
      Ropt = ', parameter_subset$MAP[parameter_subset$trait == "growth_rate" & parameter_subset$parameter == "Ropt"], '))
}')))

  viability_rate <<- eval(parse(text = paste0('function(x) {
  return(viability_rate_hourly(x,
                               Tmin = ', parameter_subset$MAP[parameter_subset$trait == "development_rate" & parameter_subset$parameter == "Tmin"], ',
                               Tmax = ', parameter_subset$MAP[parameter_subset$trait == "development_rate" & parameter_subset$parameter == "Tmax"], ',
                               Topt = ', parameter_subset$MAP[parameter_subset$trait == "development_rate" & parameter_subset$parameter == "Topt"], ',
                               Ropt = ', parameter_subset$MAP[parameter_subset$trait == "development_rate" & parameter_subset$parameter == "Ropt"], ',
                               a = ', parameter_subset$MAP[parameter_subset$trait == "viability" & parameter_subset$parameter == "a"], ',
                               h = ', parameter_subset$MAP[parameter_subset$trait == "viability" & parameter_subset$parameter == "h"], ',
                               k = ', parameter_subset$MAP[parameter_subset$trait == "viability" & parameter_subset$parameter == "k"], '))
}')))
}
