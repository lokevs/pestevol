#### Load packages, functions, and data ----
source("packages_functions_data.R")

#### Define rate functions ----

invlogit <- function(x) (1 / (1 + exp(-x)))
fecmax_scaling_anc <- function(x) {3.35942108123634 + 0.357248526159367 * log(x)}
fecmax_scaling_cold <- function(x) {3.34826094488102 + 0.357248526159367 * log(x)}
fecmax_scaling_hot <- function(x) {3.2352203671293 + 0.357248526159367 * log(x)}

LRF_hourly <- function(temp, Tmin, Tmax, Topt, Ropt) {
  valid_idx <- temp >= Tmin & temp <= Tmax
  res <- numeric(length(temp))
  if (!any(valid_idx)) return(res)
  temp_valid <- temp[valid_idx] + 273.15
  Tmin_k <- Tmin + 273.15; Tmax_k <- Tmax + 273.15; Topt_k <- Topt + 273.15
  numerator <- (Ropt / 24) * (temp_valid - Tmax_k) * (temp_valid - Tmin_k)^2
  denominator <- (Topt_k - Tmin_k) * ((Topt_k - Tmin_k) * (temp_valid - Topt_k) - (Topt_k - Tmax_k) * (Topt_k + Tmin_k - 2 * temp_valid))
  res[valid_idx] <- pmax(0, ifelse(abs(denominator) < 1e-9, 0, numerator / denominator))
  return(res)
}

viability_rate_hourly <- function(temp, Tmin, Tmax, Topt, Ropt, a, h, k) {
  dev_rate <- LRF_hourly(pmin(pmax(temp, 17), 37), Tmin, Tmax, Topt, Ropt)
  prop_viable <- 1 - invlogit(a * (temp - h)^2 + k)
  return(prop_viable^dev_rate)
}

fec_function <- function(temp, Topt, fecmax, breadth, fatness) {
  20 ^ (-abs((temp - Topt) / (breadth / 2))^(fatness + 1)) * fecmax
}

oviposition_rate_anc <- function(temp, LRS) { exp(-6.1004015 + 0.06742549 * temp) * LRS }
oviposition_rate_sel <- function(temp, LRS) { exp(-5.5043781 + 0.04348192 * temp) * LRS }

## Regime x origin rate function lookup
get_rate_functions <- function(regime, origin) {

  ## Development
  if(regime=="anc" & origin=="bra") dev <- function(x) LRF_hourly(x, 10.725, 42.989, 34.091, 0.0577)
  else if(regime=="anc" & origin=="ca") dev <- function(x) LRF_hourly(x, 9.315, 42.303, 34.698, 0.0601)
  else if(regime=="anc" & origin=="yem") dev <- function(x) LRF_hourly(x, 11.145, 42.506, 34.093, 0.0606)
  else if(regime=="cold" & origin=="bra") dev <- function(x) LRF_hourly(x, 10.963, 44.593, 34.544, 0.0576)
  else if(regime=="cold" & origin=="ca") dev <- function(x) LRF_hourly(x, 10.077, 41.701, 34.050, 0.0597)
  else if(regime=="cold" & origin=="yem") dev <- function(x) LRF_hourly(x, 11.159, 43.297, 34.057, 0.0586)
  else if(regime=="hot" & origin=="bra") dev <- function(x) LRF_hourly(x, 11.276, 42.794, 34.029, 0.0560)
  else if(regime=="hot" & origin=="ca") dev <- function(x) LRF_hourly(x, 12.415, 46.064, 34.388, 0.0548)
  else if(regime=="hot" & origin=="yem") dev <- function(x) LRF_hourly(x, 11.879, 45.476, 34.401, 0.0558)

  ## Viability
  if(regime=="anc" & origin=="bra") via <- function(x) viability_rate_hourly(x, 10.725, 42.989, 34.091, 0.0577, 0.0294, 28.130, -4.654)
  else if(regime=="anc" & origin=="ca") via <- function(x) viability_rate_hourly(x, 9.315, 42.303, 34.698, 0.0601, 0.0294, 28.130, -4.654)
  else if(regime=="anc" & origin=="yem") via <- function(x) viability_rate_hourly(x, 11.145, 42.506, 34.093, 0.0606, 0.0294, 28.130, -4.654)
  else if(regime=="cold" & origin=="bra") via <- function(x) viability_rate_hourly(x, 10.963, 44.593, 34.544, 0.0576, 0.0226, 26.798, -3.918)
  else if(regime=="cold" & origin=="ca") via <- function(x) viability_rate_hourly(x, 10.077, 41.701, 34.050, 0.0597, 0.0226, 26.798, -3.918)
  else if(regime=="cold" & origin=="yem") via <- function(x) viability_rate_hourly(x, 11.159, 43.297, 34.057, 0.0586, 0.0226, 26.798, -3.918)
  else if(regime=="hot" & origin=="bra") via <- function(x) viability_rate_hourly(x, 11.276, 42.794, 34.029, 0.0560, 0.0288, 27.427, -5.482)
  else if(regime=="hot" & origin=="ca") via <- function(x) viability_rate_hourly(x, 12.415, 46.064, 34.388, 0.0548, 0.0288, 27.427, -5.482)
  else if(regime=="hot" & origin=="yem") via <- function(x) viability_rate_hourly(x, 11.879, 45.476, 34.401, 0.0558, 0.0288, 27.427, -5.482)

  ## Growth
  if(regime=="anc" & origin=="bra") gro <- function(x) LRF_hourly(x, 11.100, 44.095, 33.005, 7.80e-05)
  else if(regime=="anc" & origin=="ca") gro <- function(x) LRF_hourly(x, 8.510, 43.001, 34.153, 8.25e-05)
  else if(regime=="anc" & origin=="yem") gro <- function(x) LRF_hourly(x, 12.631, 44.012, 31.520, 7.42e-05)
  else if(regime=="cold" & origin=="bra") gro <- function(x) LRF_hourly(x, 11.198, 45.704, 33.508, 7.84e-05)
  else if(regime=="cold" & origin=="ca") gro <- function(x) LRF_hourly(x, 10.750, 42.311, 33.607, 8.91e-05)
  else if(regime=="cold" & origin=="yem") gro <- function(x) LRF_hourly(x, 11.692, 44.823, 33.196, 8.66e-05)
  else if(regime=="hot" & origin=="bra") gro <- function(x) LRF_hourly(x, 10.455, 44.547, 32.466, 9.94e-05)
  else if(regime=="hot" & origin=="ca") gro <- function(x) LRF_hourly(x, 13.958, 46.456, 33.315, 9.52e-05)
  else if(regime=="hot" & origin=="yem") gro <- function(x) LRF_hourly(x, 12.992, 46.088, 32.645, 9.40e-05)

  ## Oviposition
  if(regime=="anc") {
    if(origin=="bra") { Topt=28.14; fecmax=70.30; breadth=21.31; fatness=4.34 }
    if(origin=="ca")  { Topt=29.01; fecmax=60.37; breadth=21.51; fatness=4.34 }
    if(origin=="yem") { Topt=27.98; fecmax=68.00; breadth=22.15; fatness=4.34 }
    ovi <- function(x) oviposition_rate_anc(x, fec_function(x, Topt, fecmax, breadth, fatness))
  } else if(regime=="cold") {
    if(origin=="bra") { Topt=27.86; fecmax=62.24; breadth=21.96; fatness=5.48 }
    if(origin=="ca")  { Topt=28.16; fecmax=67.26; breadth=24.10; fatness=5.48 }
    if(origin=="yem") { Topt=28.23; fecmax=65.83; breadth=22.29; fatness=5.48 }
    ovi <- function(x) oviposition_rate_sel(x, fec_function(x, Topt, fecmax, breadth, fatness))
  } else if(regime=="hot") {
    if(origin=="bra") { Topt=29.03; fecmax=81.82; breadth=20.86; fatness=5.48 }
    if(origin=="ca")  { Topt=28.93; fecmax=81.25; breadth=22.29; fatness=5.48 }
    if(origin=="yem") { Topt=28.98; fecmax=79.78; breadth=23.92; fatness=5.48 }
    ovi <- function(x) oviposition_rate_sel(x, fec_function(x, Topt, fecmax, breadth, fatness))
  }

  ## Fecundity scaling and max
  if(regime=="anc") { scale <- fecmax_scaling_anc; fmax_val <- if(origin=="bra") 70.30 else if(origin=="ca") 60.37 else 68.00 }
  if(regime=="cold") { scale <- fecmax_scaling_cold; fmax_val <- if(origin=="bra") 62.24 else if(origin=="ca") 67.26 else 65.83 }
  if(regime=="hot") { scale <- fecmax_scaling_hot; fmax_val <- if(origin=="bra") 81.82 else if(origin=="ca") 81.25 else 79.78 }

  return(list(dev=dev, via=via, gro=gro, ovi=ovi, scale=scale, fmax=fmax_val))
}

#### Setup data and variables ----

if (!exists("simulated_time_series")) simulated_time_series <- read.delim("../output/simulated_time_series.txt")
simulated_time_series$temp <- round(simulated_time_series$temp * 2)/2

## Mortality settings
lifespan_long <- 20
lifespan_med <- 10
lifespan_short <- 5

## Adult hourly survival rates
rate_adult_high <- exp(-(1 / (lifespan_long * 24)))
rate_adult_med  <- exp(-(1 / (lifespan_med * 24)))
rate_adult_low  <- exp(-(1 / (lifespan_short * 24)))

## Juvenile hourly survival rates (2x adult lifespan)
rate_juv_high <- exp(-(1 / (lifespan_long * 2 * 24)))
rate_juv_med  <- exp(-(1 / (lifespan_med * 2 * 24)))
rate_juv_low  <- exp(-(1 / (lifespan_short * 2 * 24)))

simulated_fitness_list <- list()
iter <- 1

#### Simulation loop ----

if (!file.exists("../output/simu_dat_long.txt")) {
for(years in c(2, 77)){
  for(wr in c(0.02, 0.04, 0.06)){

    if(years == 2 & wr != 0.02) next # skip warming rates for baseline year

    current_window <- subset(simulated_time_series, year <= years + 2 & year >= years - 1 & warming.rate == wr)

    for(si in unique(current_window$site)){

      simu_site <- subset(current_window, site == si)
      simu_site <- simu_site[order(simu_site$year, simu_site$day, simu_site$hour), ]
      temp_vec <- simu_site$temp
      n_temps <- length(temp_vec)

      for(regime in c("anc", "hot", "cold")){
        for(orig in c("yem", "ca", "bra")){
          for(tr in c("all", "fecundity", "viability", "development", "growth")){

            ## Trait isolation logic
            if(tr != "all" & regime == "anc") next # only separate traits for evolved regimes

            r_fec <- regime
            r_via <- regime
            r_dev <- regime
            r_gro <- regime

            # Isolation: target trait stays evolved, others forced to ancestral
            if(tr == "fecundity"){ r_via <- "anc"; r_dev <- "anc"; r_gro <- "anc" }
            if(tr == "viability"){ r_fec <- "anc"; r_dev <- "anc"; r_gro <- "anc" }
            if(tr == "development"){ r_fec <- "anc"; r_via <- "anc"; r_gro <- "anc" }
            if(tr == "growth"){ r_fec <- "anc"; r_via <- "anc"; r_dev <- "anc" }

            ## Get rate functions (mix and match by trait isolation)
            funcs_dev <- get_rate_functions(r_dev, orig)
            funcs_via <- get_rate_functions(r_via, orig)
            funcs_gro <- get_rate_functions(r_gro, orig)
            funcs_fec <- get_rate_functions(r_fec, orig)

            dev_rate_fun <- funcs_dev$dev
            via_rate_fun <- funcs_via$via
            growth_rate_fun <- funcs_gro$gro
            ovi_rate_fun <- funcs_fec$ovi
            fec_scale_fun <- funcs_fec$scale
            fecmax_val <- funcs_fec$fmax

            ## Initialize output dataframe
            simulated_fitness <- as.data.frame(matrix(nrow = 365, ncol = 23))
            colnames(simulated_fitness) <- c("year", "site", "warming.rate", "trait", "selection.regime", "origin", "day",
                                             "dev.time", "offspring.mass", "starting.mass", "viability",
                                             "offspring.high.mort", "offspring.med.mort", "offspring.low.mort",
                                             "total.offspring.mass.high.mort", "total.offspring.mass.med.mort", "total.offspring.mass.low.mort",
                                             "offspring.growth.rate.high.mort", "offspring.growth.rate.med.mort", "offspring.growth.rate.low.mort",
                                             "fitness.rate.high.mort", "fitness.rate.med.mort", "fitness.rate.low.mort")

            simulated_fitness$day <- 1:365
            simulated_fitness$year <- years + 2023
            simulated_fitness$warming.rate <- if(years == 2) 0 else wr
            simulated_fitness$trait <- tr
            simulated_fitness$selection.regime <- regime
            simulated_fitness$origin <- orig
            simulated_fitness$site <- si

            ## Vectorized rate calculations
            all_dev <- dev_rate_fun(temp_vec)
            all_gro <- growth_rate_fun(temp_vec)
            all_via <- via_rate_fun(temp_vec)
            all_ovi <- ovi_rate_fun(temp_vec)
            all_via_log <- log(all_via + 1e-15)

            ## Cumulative sums
            C_dev <- c(0, cumsum(all_dev))
            C_gro <- c(0, cumsum(all_gro))
            C_via <- c(0, cumsum(all_via_log))
            C_ovi <- c(0, cumsum(all_ovi))

            ## Development time (forward)
            main_cohorts <- 1:(365 + 50) # Pad to handle late-year oviposition
            main_start_h <- (364 + main_cohorts) * 24 + 1

            dev_targets <- 1 + C_dev[main_start_h]
            dev_end_h <- findInterval(dev_targets, C_dev)

            dev_h <- dev_end_h - main_start_h + 1
            dev_h[dev_end_h >= n_temps] <- NA
            incomplete <- is.na(dev_h)

            dev_h[incomplete] <- 365 * 24
            dev_h[!incomplete & dev_h > 364*24] <- 365 * 24

            simulated_fitness$dev.time <- dev_h[1:365] / 24

            ## Offspring mass
            mass_vec <- C_gro[dev_end_h + 1] - C_gro[main_start_h]
            mass_vec[incomplete] <- 0
            simulated_fitness$offspring.mass <- mass_vec[1:365]

            ## Viability
            via_vec <- exp(C_via[dev_end_h + 1] - C_via[main_start_h])
            via_vec[is.na(via_vec)] <- 0
            simulated_fitness$viability <- via_vec[1:365]

            ## Starting mass (backward)
            rev_dev <- rev(all_dev); rev_gro <- rev(all_gro)
            C_rev_dev <- c(0, cumsum(rev_dev)); C_rev_gro <- c(0, cumsum(rev_gro))
            rev_start <- n_temps - main_start_h + 1
            rev_targ <- 1 + C_rev_dev[rev_start]
            rev_end <- findInterval(rev_targ, C_rev_dev)

            start_mass <- C_rev_gro[rev_end + 1] - C_rev_gro[rev_start]
            start_mass[is.na(start_mass)] <- 0
            simulated_fitness$starting.mass <- start_mass[1:365]

            ## Daily loop for mortality-specific metrics
            for(i in 1:365) {
              if (simulated_fitness$starting.mass[i] <= 0) next
              
              start_h <- main_start_h[i]
              fec_scaled <- fecmax_val * fec_scale_fun(simulated_fitness$starting.mass[i])
              
              # High mortality (short lifespan)
              max_cnt_high <- min(which(cumsum(all_ovi[start_h:n_temps]) >= fec_scaled)[1],
                                  which(rate_adult_low ^ (1:10000) <= 0.01)[1], na.rm = T)
              if(is.na(max_cnt_high)) max_cnt_high <- 1
              ovi_lst_high <- all_ovi[start_h:(start_h + max_cnt_high - 1)] * rate_adult_low^(1:max_cnt_high)
              day_off_high <- colSums(matrix(c(ovi_lst_high, rep(0, 24 - length(ovi_lst_high) %% 24)), 24))
              
              idx_high <- i:(i + length(day_off_high) - 1)
              gen_time_high <- (1:length(day_off_high)) + (dev_h[idx_high] / 24)
              P_juv_high <- rate_juv_low ^ dev_h[idx_high]
              
              simulated_fitness$offspring.high.mort[i] <- sum(day_off_high * 0.5 * P_juv_high * via_vec[idx_high], na.rm = TRUE)
              simulated_fitness$total.offspring.mass.high.mort[i] <- sum(day_off_high * mass_vec[idx_high] * P_juv_high * via_vec[idx_high], na.rm = TRUE)
              simulated_fitness$offspring.growth.rate.high.mort[i] <- sum(day_off_high * mass_vec[idx_high] / gen_time_high * P_juv_high * via_vec[idx_high], na.rm = TRUE)
              simulated_fitness$fitness.rate.high.mort[i] <- sum(day_off_high * 0.5 / gen_time_high * P_juv_high * via_vec[idx_high], na.rm = TRUE)
              
              # Medium mortality (intermediate lifespan)
              max_cnt_med <- min(which(cumsum(all_ovi[start_h:n_temps]) >= fec_scaled)[1],
                                 which(rate_adult_med ^ (1:10000) <= 0.01)[1], na.rm = T)
              if(is.na(max_cnt_med)) max_cnt_med <- 1
              ovi_lst_med <- all_ovi[start_h:(start_h + max_cnt_med - 1)] * rate_adult_med^(1:max_cnt_med)
              day_off_med <- colSums(matrix(c(ovi_lst_med, rep(0, 24 - length(ovi_lst_med) %% 24)), 24))
              
              idx_med <- i:(i + length(day_off_med) - 1)
              gen_time_med <- (1:length(day_off_med)) + (dev_h[idx_med] / 24)
              P_juv_med <- rate_juv_med ^ dev_h[idx_med]
              
              simulated_fitness$offspring.med.mort[i] <- sum(day_off_med * 0.5 * P_juv_med * via_vec[idx_med], na.rm = TRUE)
              simulated_fitness$total.offspring.mass.med.mort[i] <- sum(day_off_med * mass_vec[idx_med] * P_juv_med * via_vec[idx_med], na.rm = TRUE)
              simulated_fitness$offspring.growth.rate.med.mort[i] <- sum(day_off_med * mass_vec[idx_med] / gen_time_med * P_juv_med * via_vec[idx_med], na.rm = TRUE)
              simulated_fitness$fitness.rate.med.mort[i] <- sum(day_off_med * 0.5 / gen_time_med * P_juv_med * via_vec[idx_med], na.rm = TRUE)
              
              # Low mortality (long lifespan)
              max_cnt_low <- min(which(cumsum(all_ovi[start_h:n_temps]) >= fec_scaled)[1],
                                 which(rate_adult_high ^ (1:10000) <= 0.01)[1], na.rm = T)
              if(is.na(max_cnt_low)) max_cnt_low <- 1
              ovi_lst_low <- all_ovi[start_h:(start_h + max_cnt_low - 1)] * rate_adult_high^(1:max_cnt_low)
              day_off_low <- colSums(matrix(c(ovi_lst_low, rep(0, 24 - length(ovi_lst_low) %% 24)), 24))
              
              idx_low <- i:(i + length(day_off_low) - 1)
              gen_time_low <- (1:length(day_off_low)) + (dev_h[idx_low] / 24)
              P_juv_low <- rate_juv_high ^ dev_h[idx_low]
              
              simulated_fitness$offspring.low.mort[i] <- sum(day_off_low * 0.5 * P_juv_low * via_vec[idx_low], na.rm = TRUE)
              simulated_fitness$total.offspring.mass.low.mort[i] <- sum(day_off_low * mass_vec[idx_low] * P_juv_low * via_vec[idx_low], na.rm = TRUE)
              simulated_fitness$offspring.growth.rate.low.mort[i] <- sum(day_off_low * mass_vec[idx_low] / gen_time_low * P_juv_low * via_vec[idx_low], na.rm = TRUE)
              simulated_fitness$fitness.rate.low.mort[i] <- sum(day_off_low * 0.5 / gen_time_low * P_juv_low * via_vec[idx_low], na.rm = TRUE)
              }

            simulated_fitness_list[[iter]] <- simulated_fitness
            print(paste("Iter:", iter, "| Year:", years, "| Regime:", regime, "| Trait:", tr))
            iter <- iter + 1
          }
        }
      }
    }
  }
}

#### Reshape and export simulation results ----

library(tidyr)
library(dplyr)
simu_dat <- do.call(rbind, simulated_fitness_list)

simu_dat_long <- simu_dat %>%
  pivot_longer(
    cols = matches("mort"),
    names_to = c(".value", "mortality"),
    names_pattern = "(.*)\\.(high|med|low)\\.mort$"
  )

write.table(simu_dat_long, "../output/simu_dat_long.txt", row.names = FALSE, quote = FALSE, sep = "\t")
} # end if !file.exists

#### Figures ----

simu_dat <- subset(read.delim("../output/simu_dat_long.txt"), mortality == "med")

## Winter vs summer fitness
pdf("../figures/raw/fig_3_winter_vs_summer_fitness.pdf", height = 1.5, width = 3, pointsize = 3); {

  simu_dat_aggregated_days <- merge(aggregate(fitness.rate ~ trait + selection.regime + day + year + warming.rate, data = simu_dat , "mean"),
                                    aggregate(offspring.growth.rate ~ trait + selection.regime + day + year + warming.rate, data = simu_dat, "mean"),
                                    by = c("trait", "selection.regime", "day","year", "warming.rate"))

  fitness_change_cold<- ifelse(subset(simu_dat_aggregated_days, selection.regime == "cold" & trait == "all" & warming.rate == 0)$fitness.rate > subset(simu_dat_aggregated_days, selection.regime == "hot" & trait == "all" & warming.rate == 0)$fitness.rate,
                               (subset(simu_dat_aggregated_days, selection.regime == "cold" & trait == "all" & warming.rate == 0)$fitness.rate / subset(simu_dat_aggregated_days, selection.regime == "hot" & trait == "all" & warming.rate == 0)$fitness.rate - 1) * 100,
                               0)

  fitness_change_warm <- ifelse(subset(simu_dat_aggregated_days, selection.regime == "cold" & trait == "all" & warming.rate == 0)$fitness.rate > subset(simu_dat_aggregated_days, selection.regime == "hot" & trait == "all" & warming.rate == 0)$fitness.rate,
                                0,
                                (subset(simu_dat_aggregated_days, selection.regime == "hot" & trait == "all" & warming.rate == 0)$fitness.rate / subset(simu_dat_aggregated_days, selection.regime == "cold" & trait == "all" & warming.rate == 0)$fitness.rate - 1 )* 100)

  days_vector <- subset(simu_dat_aggregated_days, selection.regime == "cold" & trait == "all" & warming.rate == 0)$day

  plot(NA, xlim = c(0, 366), ylim = c(0, 100), ylab = "", yaxt = "n", bty = "l", xaxt = "n", xlab = "")
  polygon(fitness_change_warm[order(days_vector)], col = alpha("orangered", 0.5), border = F)
  polygon(c(0, fitness_change_cold[order(days_vector)][-1], 0), border = F, col = alpha("dodgerblue", 0.5))
  abline(h = 0, lwd = 2, col = "white")
  axis(2, at = c(0, 50, 100), labels = c("0%", "50%", "100%"))
  axis(1, at = c(1, 365), labels = c("Jan 1", "Dec 31"))
  abline(v = 87) # March 27
  abline(v = 276) # Oct 3

}; dev.off()

## Three example sites
pdf("../figures/raw/fig_3_three_example_sites.pdf", height = 1.5, width = 1.5, pointsize = 3); {
  par(mfrow = c(1,1))
  plot(NA, xlim = c(0, 366), ylim = c(0, 0.75), bty = "l", main = "McArthur", xaxt = "n", yaxt = "n", ylab = "", xlab = "", bty = "n")
  axis(1, c(1, 365), labels = c("Jan 1", "Dec 31"))
  legend("topleft", bty = "n", paste0("+", round(((sum(subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="McArthur" & warming.rate == 0.04)$fitness.rate)/
                                                     sum(subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="McArthur" & warming.rate == 0)$fitness.rate))-1) * 100), "%"))
  polygon(c(0, aggregate(fitness.rate ~ day, subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="McArthur" & warming.rate == 0.04), "mean")$fitness.rate, 0), col = alpha("lightblue", 0.75), border = F)
  polygon(c(0, aggregate(fitness.rate ~ day, subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="McArthur" & warming.rate == 0), "mean")$fitness.rate, 0), col = alpha("blue4", 0.75), border = F)
  abline(h = 0, col = "white")
  plot(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.04 & simulated_time_series$site == "McArthur" & simulated_time_series$year == 77],
       xaxt = "n", yaxt = "n", ylab = "", xlab = "", bty = "n", type = "l", bty = "l", col = alpha("lightblue", 0.75), lwd = 0.5, ylim = c(-10, 50))
  legend("topleft", bty = "n", paste0(round(mean(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.02 & simulated_time_series$site == "McArthur" & simulated_time_series$year == 1])),
                                      " - ", round(mean(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.04 & simulated_time_series$site == "McArthur" & simulated_time_series$year == 77])),
                                      "°C mean"))
  lines(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.02 & simulated_time_series$site == "McArthur" & simulated_time_series$year == 1],
        ylab = "", xlab = "", type = "l", bty = "l", col = alpha("blue4",0.75), lwd = 0.5)
  axis(1, c(1, 365 * 24), labels = c("Jan 1", "Dec 31"))
  axis(2, c(0, 40), las = 2)

  plot(NA, xlim = c(0, 366), ylim = c(0, 0.75), bty = "l", main = "FivePoints", xaxt = "n", yaxt = "n", ylab = "", xlab = "", bty = "n")
  axis(1, c(1, 365), labels = c("Jan 1", "Dec 31"))
  legend("topleft", bty = "n", paste0("+", round(((sum(subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="FivePoints" & warming.rate == 0.04)$fitness.rate)/
                                                     sum(subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="FivePoints" & warming.rate == 0)$fitness.rate))-1) * 100), "%"))
  polygon(c(0, aggregate(fitness.rate ~ day, subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="FivePoints" & warming.rate == 0.04), "mean")$fitness.rate, 0), col = alpha("yellow2", 0.75), border = F)
  polygon(c(0, aggregate(fitness.rate ~ day, subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="FivePoints" & warming.rate == 0), "mean")$fitness.rate, 0), col = alpha("orange", 0.75), border = F)
  abline(h = 0, col = "white")
  plot(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.04 & simulated_time_series$site == "FivePoints" & simulated_time_series$year == 77],
       ylab = "", xlab = "", yaxt = "n", xaxt = "n", type = "l", bty = "l", col = alpha("yellow2", 0.75), lwd = 0.5, ylim = c(-10, 50))
  legend("topleft", bty = "n", paste0(round(mean(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.02 & simulated_time_series$site == "FivePoints" & simulated_time_series$year == 1])),
                                      " - ", round(mean(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.04 & simulated_time_series$site == "FivePoints" & simulated_time_series$year == 77])),
                                      "°C mean"))
  lines(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.02 & simulated_time_series$site == "FivePoints" & simulated_time_series$year == 1],
        type = "l", col = alpha("orange", 0.75), lwd = 0.5)
  axis(1, c(1, 365 * 24), labels = c("Jan 1", "Dec 31"))
  axis(2, c(0, 40), las = 2)

  plot(NA, xlim = c(0, 366), ylim = c(0, 0.75), bty = "l", main = "Seeley", xaxt = "n", yaxt = "n", ylab = "", xlab = "", bty = "n")
  axis(1, c(1, 365), labels = c("Jan 1", "Dec 31"))
  legend("topleft", bty = "n", paste0("+", round(((sum(subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="Seeley" & warming.rate == 0.04)$fitness.rate)/
                                                     sum(subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="Seeley" & warming.rate == 0)$fitness.rate))-1) * 100), "%"))
  polygon(c(0, aggregate(fitness.rate ~ day, subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="Seeley" & warming.rate == 0.04), "mean")$fitness.rate, 0), col = alpha("pink", 0.75), border = F)
  polygon(c(0, aggregate(fitness.rate ~ day, subset(simu_dat, selection.regime == "anc" & trait == "all" & site =="Seeley" & warming.rate == 0), "mean")$fitness.rate, 0), col = alpha("red2", 0.75), border = F)
  abline(h = 0, col = "white")
  plot(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.04 & simulated_time_series$site == "Seeley" & simulated_time_series$year == 77],
       ylab = "", xlab = "", yaxt = "n", xaxt = "n", type = "l", bty = "l", col = alpha("pink", 0.75), lwd = 0.5, ylim = c(-10, 50))
  legend("topleft", bty = "n", paste0(round(mean(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.02 & simulated_time_series$site == "Seeley" & simulated_time_series$year == 1])),
                                      " - ", round(mean(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.04 & simulated_time_series$site == "Seeley" & simulated_time_series$year == 77])),
                                      "°C mean"))
  lines(simulated_time_series$temperature[simulated_time_series$warming.rate == 0.02 & simulated_time_series$site == "Seeley" & simulated_time_series$year == 1],
        type = "l", col = alpha("red", 0.75), lwd = 0.5)
  axis(1, c(1, 365 * 24), labels = c("Jan 1", "Dec 31"))
  axis(2, c(0, 40), las = 2)

}; dev.off()

## Seasonal absolute fitness differences (baseline and warmed)
pdf("../figures/raw/fig_3_seasonal_absolute_fitness_differences.pdf", height = 1.5, width = 1.5, pointsize = 3); {

  simu_dat$batch <- ceiling(simu_dat$day / 5)
  aggregated_fitness <- aggregate(fitness.rate ~ day + trait + origin + selection.regime + site + warming.rate,
                                  data = simu_dat,
                                  FUN = mean)

  evolved_lines_df <- aggregated_fitness[aggregated_fitness$selection.regime != "anc", ]
  ancestral_lines_df <- aggregated_fitness[aggregated_fitness$selection.regime == "anc", ]
  ancestral_baseline <- ancestral_lines_df[, c("day", "trait", "site", "warming.rate", "origin", "fitness.rate")]
  names(ancestral_baseline)[names(ancestral_baseline) == "fitness.rate"] <- "anc_fitness_rate"
  aggregated_fitness_evolved <- merge(evolved_lines_df, ancestral_baseline,
                                      by = c("day", "trait", "origin", "site", "warming.rate"),
                                      all.x = TRUE)
  aggregated_fitness_evolved$fitness.change <- aggregated_fitness_evolved$fitness.rate - aggregated_fitness_evolved$anc_fitness_rate

  # Baseline (warming.rate == 0)
  plot(NA, ylim = c(-0.1, 0.16), xlim = c(0, 366), xaxt = "n", yaxt = "n", ylab = "rmax difference", xlab = "", bty = "l")
  abline(h = 0, lty = 3)
  axis(1, at = c(0, 365), c("Jan 1", "Dec 31"))
  for(regime in c("cold", "hot")){
    for(ori in c("bra", "yem", "ca")){
      if(regime == "cold")(colour <- "dodgerblue")
      else(colour <- "orangered")

      line <- aggregate(fitness.change ~ day, subset(aggregated_fitness_evolved, selection.regime == regime & origin == ori & warming.rate == 0), "mean")
      interval <- aggregate(fitness.change ~ day, subset(aggregated_fitness_evolved, selection.regime == regime & origin == ori & warming.rate == 0), "sd")

      polygon(y = c(line[,2] + interval[,2], rev(line[,2] - interval[,2])),
              x = c(line[,1], rev(line[,1])),
              border = FALSE,
              col = alpha(colour, 0.4))
    }
  }
  for(regime in c("cold", "hot")){
    for(ori in c("bra", "yem", "ca")){
      if(regime == "cold")(colour <- "dodgerblue")
      else(colour <- "orangered")
      if(ori == "bra")(linetype <- 1)
      if(ori == "yem")(linetype <-  2)
      if(ori == "ca")(linetype <-  3)

      line <- aggregate(fitness.change ~ day, subset(aggregated_fitness_evolved, selection.regime == regime & origin == ori & warming.rate == 0), "mean")
      lines(line, col = colour, lwd = 1, lty = linetype)
    }
  }

  # Warmed (warming.rate == 0.04)
  plot(NA, ylim = c(-0.1, 0.16), xlim = c(0, 366), xaxt = "n", yaxt = "n", ylab = "rmax difference", xlab = "", bty = "l")
  abline(h = 0, lty = 3)
  axis(1, at = c(0, 365), c("Jan 1", "Dec 31"))
  for(regime in c("cold", "hot")){
    for(ori in c("bra", "yem", "ca")){
      if(regime == "cold")(colour <- "dodgerblue")
      else(colour <- "orangered")

      line <- aggregate(fitness.change ~ day, subset(aggregated_fitness_evolved, selection.regime == regime & origin == ori & warming.rate == 0.04), "mean")
      interval <- aggregate(fitness.change ~ day, subset(aggregated_fitness_evolved, selection.regime == regime & origin == ori & warming.rate == 0.04), "sd")

      polygon(y = c(line[,2] + interval[,2], rev(line[,2] - interval[,2])),
              x = c(line[,1], rev(line[,1])),
              border = FALSE,
              col = alpha(colour, 0.4))
    }
  }
  for(regime in c("cold", "hot")){
    for(ori in c("bra", "yem", "ca")){
      if(regime == "cold")(colour <- "dodgerblue")
      else(colour <- "orangered")
      if(ori == "bra")(linetype <- 1)
      if(ori == "yem")(linetype <-  2)
      if(ori == "ca")(linetype <-  3)

      line <- aggregate(fitness.change ~ day, subset(aggregated_fitness_evolved, selection.regime == regime & origin == ori & warming.rate == 0.04), "mean")
      lines(line, col = colour, lwd = 1, lty = linetype)
    }
  }

}; dev.off()

## Seasonal absolute crop damage differences (baseline and warmed)
pdf("../figures/raw/fig_3_seasonal_absolute_crop_damage_differences.pdf", height = 1.5, width = 1.5, pointsize = 3); {

  simu_dat$batch <- ceiling(simu_dat$day / 5)
  aggregated_crop_damage <- aggregate(offspring.growth.rate ~ day + trait + origin + selection.regime + site + warming.rate,
                                      data = simu_dat,
                                      FUN = mean)

  evolved_lines_df <- aggregated_crop_damage[aggregated_crop_damage$selection.regime != "anc", ]
  ancestral_lines_df <- aggregated_crop_damage[aggregated_crop_damage$selection.regime == "anc", ]
  ancestral_baseline <- ancestral_lines_df[, c("day", "trait", "site", "warming.rate", "origin", "offspring.growth.rate")]
  names(ancestral_baseline)[names(ancestral_baseline) == "offspring.growth.rate"] <- "anc_crop_damage_rate"
  aggregated_crop_damage_evolved <- merge(evolved_lines_df, ancestral_baseline,
                                          by = c("day", "trait", "origin", "site", "warming.rate"),
                                          all.x = TRUE)
  aggregated_crop_damage_evolved$crop.damage.change <- aggregated_crop_damage_evolved$offspring.growth.rate - aggregated_crop_damage_evolved$anc_crop_damage_rate

  # Baseline (warming.rate == 0)
  plot(NA, ylim = c(-0.0001, 0.0005), xlim = c(0, 366), xaxt = "n", yaxt = "n", ylab = "rmax difference", xlab = "", bty = "l")
  abline(h = 0, lty = 3)
  axis(1, at = c(0, 365), c("Jan 1", "Dec 31"))
  for(regime in c("cold", "hot")){
    for(ori in c("bra", "yem", "ca")){
      if(regime == "cold")(colour <- "dodgerblue")
      else(colour <- "orangered")

      line <- aggregate(crop.damage.change ~ day, subset(aggregated_crop_damage_evolved, selection.regime == regime & origin == ori & warming.rate == 0), "mean")
      interval <- aggregate(crop.damage.change ~ day, subset(aggregated_crop_damage_evolved, selection.regime == regime & origin == ori & warming.rate == 0), "sd")

      polygon(y = c(line[,2] + interval[,2], rev(line[,2] - interval[,2])),
              x = c(line[,1], rev(line[,1])),
              border = FALSE,
              col = alpha(colour, 0.4))
    }
  }
  for(regime in c("cold", "hot")){
    for(ori in c("bra", "yem", "ca")){
      if(regime == "cold")(colour <- "dodgerblue")
      else(colour <- "orangered")
      if(ori == "bra")(linetype <- 1)
      if(ori == "yem")(linetype <-  2)
      if(ori == "ca")(linetype <-  3)

      line <- aggregate(crop.damage.change ~ day, subset(aggregated_crop_damage_evolved, selection.regime == regime & origin == ori & warming.rate == 0), "mean")
      lines(line, col = colour, lwd = 1, lty = linetype)
    }
  }

  # Warmed (warming.rate == 0.04)
  plot(NA, ylim = c(-0.0001, 0.0005), xlim = c(0, 366), xaxt = "n", yaxt = "n", ylab = "rmax difference", xlab = "", bty = "l")
  abline(h = 0, lty = 3)
  axis(1, at = c(0, 365), c("Jan 1", "Dec 31"))
  for(regime in c("cold", "hot")){
    for(ori in c("bra", "yem", "ca")){
      if(regime == "cold")(colour <- "dodgerblue")
      else(colour <- "orangered")

      line <- aggregate(crop.damage.change ~ day, subset(aggregated_crop_damage_evolved, selection.regime == regime & origin == ori & warming.rate == 0.04), "mean")
      interval <- aggregate(crop.damage.change ~ day, subset(aggregated_crop_damage_evolved, selection.regime == regime & origin == ori & warming.rate == 0.04), "sd")

      polygon(y = c(line[,2] + interval[,2], rev(line[,2] - interval[,2])),
              x = c(line[,1], rev(line[,1])),
              border = FALSE,
              col = alpha(colour, 0.4))
    }
  }
  for(regime in c("cold", "hot")){
    for(ori in c("bra", "yem", "ca")){
      if(regime == "cold")(colour <- "dodgerblue")
      else(colour <- "orangered")
      if(ori == "bra")(linetype <- 1)
      if(ori == "yem")(linetype <-  2)
      if(ori == "ca")(linetype <-  3)

      line <- aggregate(crop.damage.change ~ day, subset(aggregated_crop_damage_evolved, selection.regime == regime & origin == ori & warming.rate == 0.04), "mean")
      lines(line, col = colour, lwd = 1, lty = linetype)
    }
  }

}; dev.off()

## Fitness fold change vs baseline fitness by site
pdf("../figures/raw/fig_3_bad_sites_gain_more.pdf", height = 1.5, width = 1.7, pointsize = 3); {
  aggregates <- aggregate(fitness.rate ~ site + warming.rate + origin + selection.regime, subset(simu_dat, trait == "all"), "mean")

  aggregates$log2.fitness.rate.fold.change <- 0
  aggregates$original.fitness.rate <- 0
  aggregates$color <- ifelse(aggregates$warming.rate == 0, "darkgreen", ifelse(aggregates$warming.rate == 0.02, "orange", ifelse(aggregates$warming.rate == 0.04, "orangered", "darkred")))
  for(i in 1:nrow(aggregates)){
    aggregates$log2.fitness.rate.fold.change[i] <- log2(aggregates$fitness.rate[i] / aggregates$fitness.rate[aggregates$warming.rate == 0 & aggregates$site == aggregates$site[i] & aggregates$origin == aggregates$origin[i] & aggregates$selection.regime == aggregates$selection.regime[i]])
    aggregates$original.fitness.rate[i] <- aggregates$fitness.rate[aggregates$warming.rate == 0 & aggregates$site == aggregates$site[i] & aggregates$origin == aggregates$origin[i] & aggregates$selection.regime == aggregates$selection.regime[i]]
  }
  aggregates <- aggregates[aggregates$warming.rate != 0,]
  aggregates <- aggregates[order(aggregates$log2.fitness.rate.fold.change, decreasing = F), ]

  # No growth season constraint
  plot(sqrt(abs(log2.fitness.rate.fold.change))*sign(log2.fitness.rate.fold.change) ~ original.fitness.rate,
       ylim = c(-sqrt(1), sqrt(4)), data = aggregates[aggregates$warming.rate != 0, ],
       main = "No growth season constraint", yaxt = "n", ylab = "sqrt(log2 fitness fold change)", cex = 1, pch = 16, col = aggregates[aggregates$warming.rate != 0, ]$color, bty = "l")
  abline(h = 0)
  axis(2, c(-sqrt(1), 0, sqrt(1), sqrt(4)), c(-1, 0, 1, 4))

  # Growth season fixed (days 90-275)
  aggregates <- aggregate(fitness.rate ~ site + warming.rate + origin + selection.regime, subset(simu_dat, trait == "all" & day > 90 & day < (365-90)), "mean")

  aggregates$log2.fitness.rate.fold.change <- 0
  aggregates$original.fitness.rate <- 0
  aggregates$color <- ifelse(aggregates$warming.rate == 0, "darkgreen", ifelse(aggregates$warming.rate == 0.02, "orange", ifelse(aggregates$warming.rate == 0.04, "orangered", "darkred")))
  for(i in 1:nrow(aggregates)){
    aggregates$log2.fitness.rate.fold.change[i] <- log2(aggregates$fitness.rate[i] / aggregates$fitness.rate[aggregates$warming.rate == 0 & aggregates$site == aggregates$site[i] & aggregates$origin == aggregates$origin[i] & aggregates$selection.regime == aggregates$selection.regime[i]])
    aggregates$original.fitness.rate[i] <- aggregates$fitness.rate[aggregates$warming.rate == 0 & aggregates$site == aggregates$site[i] & aggregates$origin == aggregates$origin[i] & aggregates$selection.regime == aggregates$selection.regime[i]]
  }
  aggregates <- aggregates[aggregates$warming.rate != 0,]
  aggregates <- aggregates[order(aggregates$log2.fitness.rate.fold.change, decreasing = F), ]
  plot(sqrt(abs(log2.fitness.rate.fold.change))*sign(log2.fitness.rate.fold.change) ~ original.fitness.rate,
       ylim = c(-sqrt(1), sqrt(4)), data = aggregates[aggregates$warming.rate != 0, ],
       main = "Growth season fixed", yaxt = "n", ylab = "sqrt(log2 fitness fold change)", cex = 1, pch = 16, col = aggregates[aggregates$warming.rate != 0, ]$color, bty = "l")
  abline(h = 0)
  axis(2, c(-sqrt(1), 0, sqrt(1), sqrt(4)), c(-1, 0, 1, 4))

}; dev.off()
