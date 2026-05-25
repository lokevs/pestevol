#### Base Functions ####

# The inverse logit function 
invlogit <- function(x) (1 / (1 + exp(-x)))

# Fecundity ~ body size scaling functions (Retained from original script)
fecmax_scaling_anc <- function(x) {3.35942108123634 + 0.357248526159367 * log(x)}
fecmax_scaling_cold <- function(x) {3.34826094488102 + 0.357248526159367 * log(x)}
fecmax_scaling_hot <- function(x) {3.2352203671293 + 0.357248526159367 * log(x)}

# LRF TPC function
LRF_hourly <- function(temp, Tmin, Tmax, Topt, Ropt) {
  valid_idx <- temp >= Tmin & temp <= Tmax
  res <- numeric(length(temp))
  if (!any(valid_idx)) return(res)
  
  temp_valid <- temp[valid_idx]
  temp_k <- temp_valid + 273.15
  Tmin_k <- Tmin + 273.15
  Tmax_k <- Tmax + 273.15
  Topt_k <- Topt + 273.15
  
  numerator <- (Ropt / 24) * (temp_k - Tmax_k) * (temp_k - Tmin_k)^2
  term1 <- (Topt_k - Tmin_k)
  term2 <- term1 * (temp_k - Topt_k)
  term3 <- (Topt_k - Tmax_k) * (Topt_k + Tmin_k - 2 * temp_k)
  denominator <- term1 * (term2 - term3)
  
  calculated_vals <- ifelse(abs(denominator) < 1e-9, 0, numerator / denominator)
  res[valid_idx] <- pmax(0, calculated_vals)
  return(res)
}

# Viability function
viab_function <- function(temp, a, h, k) {
  a * (temp - h)^2 + k
}

# Viability rate function
viability_rate_hourly <- function(temp, Tmin, Tmax, Topt, Ropt, a, h, k) {
  dev_temp <- pmin(pmax(temp, 17), 37)
  dev_rate <- LRF_hourly(dev_temp, Tmin, Tmax, Topt, Ropt)
  prop_viable <- 1 - invlogit(viab_function(temp, a, h, k))
  return(prop_viable^dev_rate)
}

# Fecundity function
fec_function <- function(temp, Topt, fecmax, breadth, fatness) {
  20 ^ (-abs((temp - Topt) / (breadth / 2))^(fatness + 1)) * fecmax
}

# Oviposition rate functions 
oviposition_rate_function_anc_hourly <- function(temp, LRS) {
  return(exp(-6.10040153051259 + 0.0674254905014758 * temp) * LRS)
}
oviposition_rate_function_selection_hourly <- function(temp, LRS) {
  return(exp(-5.50437809675088 + 0.0434819157913075 * temp) * LRS)
}

#### Global Hourly Rate Functions ####

## Development rates (Values from parameter_estimates_dev_rate_global.txt)
anc_dev <- function(x) {
  return(LRF_hourly(x,
                    Tmin = 10.921086,
                    Tmax = 43.077166,
                    Topt = 34.212807,
                    Ropt = 0.0588373))
}

cold_dev <- function(x) {
  return(LRF_hourly(x,
                    Tmin = 10.799862,
                    Tmax = 43.320158,
                    Topt = 34.166528,
                    Ropt = 0.0579752))
}

hot_dev <- function(x) {
  return(LRF_hourly(x,
                    Tmin = 11.820465,
                    Tmax = 44.276211,
                    Topt = 34.179006,
                    Ropt = 0.0555860))
}

## Viability rates (Values from parameter_estimates_viab_fec_global.txt + dev params)
anc_via <- function(x) {
  return(viability_rate_hourly(x,
                               Tmin = 10.921086,
                               Tmax = 43.077166,
                               Topt = 34.212807,
                               Ropt = 0.0588373,
                               a = 0.0306695,
                               h = 28.740311, 
                               k = -5.171095))
}

cold_via <- function(x) {
  return(viability_rate_hourly(x,
                               Tmin = 10.799862,
                               Tmax = 43.320158,
                               Topt = 34.166528,
                               Ropt = 0.0579752,
                               a = 0.0254042,
                               h = 26.805169, 
                               k = -4.756591))
}

hot_via <- function(x) {
  return(viability_rate_hourly(x,
                               Tmin = 11.820465,
                               Tmax = 44.276211,
                               Topt = 34.179006,
                               Ropt = 0.0555860,
                               a = 0.0335644,
                               h = 28.979208, 
                               k = -5.759574))
}

## Growth rates (Values from parameter_estimates_growth_rate_global.txt)
anc_growth <- function(x) {
  return(LRF_hourly(x,
                    Tmin = 11.471628,
                    Tmax = 44.010579,
                    Topt = 33.143588,
                    Ropt = 7.7744e-05))
}

cold_growth <- function(x) {
  return(LRF_hourly(x,
                    Tmin = 10.934165,
                    Tmax = 44.213506,
                    Topt = 33.686478,
                    Ropt = 8.4699e-05))
}

hot_growth <- function(x) {
  return(LRF_hourly(x,
                    Tmin = 13.160604,
                    Tmax = 44.915729,
                    Topt = 32.751166,
                    Ropt = 0.00010061))
}

## Oviposition rates (Values from parameter_estimates_viab_fec_global.txt)
anc_ovi <- function(x) {
  LRS <- fec_function(x,
                      Topt = 28.401171,
                      fecmax = 67.349470,
                      breadth = 21.849511,
                      fatness = 3.783288)
  return(oviposition_rate_function_anc_hourly(x, LRS))
}

cold_ovi <- function(x) {
  LRS <- fec_function(x,
                      Topt = 28.328872,
                      fecmax = 66.979230,
                      breadth = 23.978002,
                      fatness = 3.324542)
  return(oviposition_rate_function_selection_hourly(x, LRS))
}

hot_ovi <- function(x) {
  LRS <- fec_function(x,
                      Topt = 28.793020,
                      fecmax = 79.435320,
                      breadth = 23.375238,
                      fatness = 4.596587)
  return(oviposition_rate_function_selection_hourly(x, LRS))
}

# Maximum fecundity globals
anc_fecmax <- 67.349470
cold_fecmax <- 66.979230
hot_fecmax <- 79.435320

#### Global Fitness/crop damage calculation function ####

calculate_fitness_crop_damage <- function(temperatures, plot = FALSE, lifespan_adult = 10, trait = "all", regime = "all") {
  
  param_sets <- list(
    list(name="anc", dev=anc_dev, growth=anc_growth, via=anc_via, ovi=anc_ovi, fec_scale=fecmax_scaling_anc, fecmax=anc_fecmax),
    list(name="cold", dev=cold_dev, growth=cold_growth, via=cold_via, ovi=cold_ovi, fec_scale=fecmax_scaling_cold, fecmax=cold_fecmax),
    list(name="hot", dev=hot_dev, growth=hot_growth, via=hot_via, ovi=hot_ovi, fec_scale=fecmax_scaling_hot, fecmax=hot_fecmax)
  )
  
  anc_params <- param_sets[[1]]
  
  if (regime != "all") {
    param_sets <- param_sets[sapply(param_sets, function(p) p$name == regime)]
  }
  
  full_results_list <- list()
  
  for (params in param_sets) {
    run_params <- anc_params
    
    # Substitute evolved traits
    if (trait == "dev") { run_params$dev <- params$dev
    } else if (trait == "via") { run_params$via <- params$via
    } else if (trait == "growth") { run_params$growth <- params$growth
    } else if (trait == "fec") {
      run_params$ovi <- params$ovi
      run_params$fec_scale <- params$fec_scale
      run_params$fecmax <- params$fecmax
    } else if (trait == "all") { run_params <- params }
    
    full_results_list[[params$name]] <- run_simulation(
      temperatures = temperatures,
      dev_rate = run_params$dev,
      growth_rate = run_params$growth,
      viability_rate = run_params$via,
      ovi_rate = run_params$ovi,
      fecmax_scaling = run_params$fec_scale,
      fecmax = run_params$fecmax,
      lifespan_adult = lifespan_adult
    )
  }
  
  # Summary stats
  summary_list <- lapply(names(full_results_list), function(name) {
    res <- full_results_list[[name]]
    
    # --- STRICT GEOMETRIC MEAN ---
    # If any lambda is 0, log(0) is -Inf, mean is -Inf, and exp(-Inf) is 0.
    # This naturally returns 0 if any cohort fails completely.
    geo_mean_lambda <- exp(mean(log(res$lambda), na.rm = TRUE))
    
    data.frame(
      regime = name, 
      trait = trait,
      mean_crop_damage_rate = mean(res$crop_damage, na.rm = TRUE),
      mean_lambda = mean(res$lambda, na.rm = TRUE),
      mean_lambda_geo = geo_mean_lambda, 
      reproductive_days = sum(res$offspring >= 2, na.rm = TRUE)
    )
  })
  results_df <- do.call(rbind, summary_list)
  
  if (plot) {
    par(mfrow = c(3, 3), mar = c(4.1, 4.1, 2.1, 1.1), oma = c(0, 0, 2, 0))
    plot_vars <- list(
      list(name = "dev_time", title = "Development Time"),
      list(name = "mass", title = "Final Mass"),
      list(name = "starting_mass", title = "Starting Mass"),
      list(name = "viability", title = "Viability"),
      list(name = "offspring_survival", title = "Offspring Survival Prob."),
      list(name = "offspring", title = "R0 (Net Reproductive Rate)"),
      list(name = "lambda", title = "Daily Lambda"), 
      list(name = "crop_damage", title = "Crop Damage Rate")
    )
    
    for (p_var in plot_vars) {
      all_data <- lapply(full_results_list, `[[`, p_var$name)
      ylim <- range(unlist(all_data), na.rm = TRUE)
      plot(NA, xlim = c(1, 365), ylim = ylim, xlab = "Day", ylab = "", main = p_var$title)
      for (name in names(full_results_list)) {
        color <- ifelse(name == "hot", "red", ifelse(name == "cold", "blue", "black"))
        lines(1:365, full_results_list[[name]][[p_var$name]], col = color)
      }
    }
    plot.new()
    legend("center", legend = c("Ancestral", "Cold", "Hot"), col = c("black", "blue", "red"), lty = 1, cex = 1.2, bty = "n")
    par(mfrow = c(1, 1))
  }
  return(results_df)
}

#### Internal Simulation Engine (Vectorized) ####
run_simulation <- function(temperatures, dev_rate, growth_rate, viability_rate, ovi_rate, fecmax_scaling, fecmax, lifespan_adult, sex_ratio = 0.5) {
  
  n_temps <- length(temperatures)
  survival_rate_adult <- exp(-(1 / (lifespan_adult * 24)))
  lifespan_juvenile <- lifespan_adult * 2
  survival_rate_juvenile <- exp(-(1 / (lifespan_juvenile * 24)))  
  
  dev_cohorts <- 1:(365 + 200)
  dev_start_hours <- (364 + dev_cohorts) * 24 + 1
  dev_start_hours <- dev_start_hours[dev_start_hours < n_temps]
  
  main_cohorts <- 1:365
  main_start_hours <- (364 + main_cohorts) * 24 + 1
  
  # Pre-calculate rates
  all_dev_rates <- dev_rate(temperatures)
  all_growth_rates <- growth_rate(temperatures)
  all_via_rates <- viability_rate(temperatures)
  all_ovi_rates <- ovi_rate(temperatures)
  all_via_rates_log <- log(all_via_rates + 1e-15)
  
  C_dev <- c(0, cumsum(all_dev_rates))
  C_growth <- c(0, cumsum(all_growth_rates))
  C_via_log <- c(0, cumsum(all_via_rates_log))
  C_ovi <- c(0, cumsum(all_ovi_rates))
  
  # Development time
  dev_targets <- 1 + C_dev[dev_start_hours]
  dev_end_hours <- findInterval(dev_targets, C_dev)
  dev_end_hours[dev_end_hours >= n_temps] <- NA
  simulated_dev_time_hours <- dev_end_hours - dev_start_hours + 1
  
  incomplete_mask <- is.na(simulated_dev_time_hours)
  simulated_dev_time_hours[incomplete_mask] <- 365 * 24
  
  long_dev_mask <- !incomplete_mask & simulated_dev_time_hours > (364 * 24)
  simulated_dev_time_hours[long_dev_mask] <- 365 * 24
  
  simulated_dev_time <- simulated_dev_time_hours / 24
  cumulative_offspring_survival <- survival_rate_juvenile^simulated_dev_time_hours
  
  # Mass
  simulated_mass <- C_growth[dev_end_hours + 1] - C_growth[dev_start_hours]
  simulated_mass[incomplete_mask] <- 0
  
  # Viability
  valid_dev_mask <- !is.na(dev_end_hours)
  log_viabilities <- rep(NA, length(dev_cohorts))
  log_viabilities[valid_dev_mask] <- C_via_log[dev_end_hours[valid_dev_mask] + 1] - C_via_log[dev_start_hours[valid_dev_mask]]
  simulated_viability <- exp(log_viabilities)
  simulated_viability[is.na(simulated_viability)] <- 0
  
  # Starting mass
  rev_dev_rates <- rev(all_dev_rates)
  rev_growth_rates <- rev(all_growth_rates)
  C_rev_dev <- c(0, cumsum(rev_dev_rates))
  C_rev_growth <- c(0, cumsum(rev_growth_rates))
  rev_start_hours <- n_temps - main_start_hours + 1
  rev_dev_targets <- 1 + C_rev_dev[rev_start_hours]
  rev_dev_end_hours <- findInterval(rev_dev_targets, C_rev_dev)
  rev_dev_end_hours[rev_dev_end_hours > n_temps] <- NA
  simulated_starting_mass <- C_rev_growth[rev_dev_end_hours + 1] - C_rev_growth[rev_start_hours]
  simulated_starting_mass[is.na(simulated_starting_mass)] <- 0
  
  # Final fitness loop
  max_life_by_survival <- ceiling(log(0.01) / log(survival_rate_adult))
  simulated_offspring <- numeric(365) # This is effectively R0
  simulated_crop_damage_rate <- numeric(365)
  simulated_lambda <- numeric(365)
  
  for (i in 1:365) {
    if (is.na(simulated_starting_mass[i]) || simulated_viability[i] == 0) next
    
    start_hour <- main_start_hours[i]
    fecmax_scaled <- fecmax * fecmax_scaling(simulated_starting_mass[i])
    ovi_target <- fecmax_scaled + C_ovi[start_hour]
    ovi_end_hour <- findInterval(ovi_target, C_ovi)
    
    lifespan_hours <- min(ovi_end_hour - start_hour + 1, max_life_by_survival, n_temps - start_hour + 1, na.rm = TRUE)
    if (is.na(lifespan_hours) || lifespan_hours <= 0) next
    
    end_hour <- start_hour + lifespan_hours - 1
    hours_lived <- 1:lifespan_hours
    ovi_list <- all_ovi_rates[start_hour:end_hour] * (survival_rate_adult^hours_lived)
    
    if(sum(ovi_list) == 0) next 
    
    num_days <- ceiling(lifespan_hours / 24)
    padded_ovi <- c(ovi_list, rep(0, num_days * 24 - lifespan_hours))
    daily_offspring <- colSums(matrix(padded_ovi, nrow = 24))
    
    offspring_cohort_indices <- i:(i + num_days - 1)
    valid_indices_mask <- offspring_cohort_indices <= length(simulated_mass) & !is.na(simulated_dev_time[offspring_cohort_indices]) & daily_offspring > 0
    if (!any(valid_indices_mask)) next
    
    valid_days <- 1:num_days
    daily_offspring_trimmed <- daily_offspring[valid_indices_mask]
    valid_offspring_cohorts <- offspring_cohort_indices[valid_indices_mask]
    
    time_to_maturity <- valid_days[valid_indices_mask] + simulated_dev_time[valid_offspring_cohorts]
    offspring_baseline_survival <- cumulative_offspring_survival[valid_offspring_cohorts]
    
    offspring_thermal_viability <- simulated_viability[valid_offspring_cohorts]
    
    # R0 calculation
    simulated_offspring[i] <- sum(daily_offspring_trimmed * sex_ratio * offspring_baseline_survival * offspring_thermal_viability, na.rm = TRUE)
    
    # Crop Damage (Biomass flux)
    simulated_crop_damage_rate[i] <- sum(daily_offspring_trimmed * simulated_mass[valid_offspring_cohorts] * offspring_baseline_survival * offspring_thermal_viability / time_to_maturity, na.rm = TRUE)
    
    # Lambda Calculation: R0 ^ (1 / GenTime)
    # GenTime is weighted average age of mother when offspring are produced
    # Age = time_to_maturity (Dev time + Adult days)
    if (simulated_offspring[i] > 0) {
      numerator_gen_time <- sum(daily_offspring_trimmed * sex_ratio * offspring_baseline_survival * offspring_thermal_viability * time_to_maturity, na.rm = TRUE)
      gen_time <- numerator_gen_time / simulated_offspring[i]
      simulated_lambda[i] <- simulated_offspring[i] ^ (1 / gen_time)
    } else {
      simulated_lambda[i] <- 0
    }
  }
  
  return(list(
    dev_time = simulated_dev_time[main_cohorts],
    mass = simulated_mass[main_cohorts],
    starting_mass = simulated_starting_mass,
    viability = simulated_viability[main_cohorts],
    offspring_survival = cumulative_offspring_survival[main_cohorts],
    offspring = simulated_offspring,
    crop_damage = simulated_crop_damage_rate,
    lambda = simulated_lambda
  ))
}

