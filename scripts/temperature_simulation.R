#### Load packages, functions, and data ----
source("packages_functions_data.R")
library(maps)
library(mapdata)  # provides the "worldHires" database used by the validation map below

#### Read and prepare data ----

dat <- read.delim("../data/climate_data_california.txt")
dat$day_of_year <- yday(dat$date.time)
dat$hour <- hour(dat$date.time)

#### Run temperature simulation ----

if (!file.exists("../output/simulated_time_series.txt")) {
  set.seed(534); simulated_time_series <- simulate_temperature(
    input_temperature = dat$temp,
    input_date_time = dat$date.time,
    input_site = dat$site,
    years = 78,
    n = 20,
    warming_rate = c(0.02, 0.04, 0.06)
  )
  write.table(simulated_time_series, "../output/simulated_time_series.txt", sep = "\t", row.names = F, quote = F)
}

if (!exists("simulated_time_series")) simulated_time_series <- read.delim("../output/simulated_time_series.txt")

#### Plot simulated yearly averages ----

averages <- aggregate(temperature ~ year + simulation + site, data = subset(simulated_time_series, warming.rate == 0.04), "mean")

pdf("../figures/raw/simulated_yearly_averages.pdf", height = 2, width = 5, pointsize = 3); {
plot(NA, ylim = c(10, 30), xlim = c(1, max(averages$year) + 16.5), xlab = "Year", ylab = "Yearly average temperature", xaxt = "n", yaxt = "n", bty = "l")
axis(side = 2, at = c(10, 20, 30), las  = 2)
axis(side = 1, at = c(1,  26, 51,  76), las  = 1, labels = c("Present", "2050", "2075", "2100"))

for(i in unique(averages$simulation)){
  lines(temperature ~ year, data = subset(averages, simulation == i), lwd = 0.5, col = "grey")
}

## Select extreme and mid regimes
mean_of_means <- aggregate(temperature ~ simulation + site, data = averages, "mean")
select_simulations <- c(mean_of_means$simulation[mean_of_means$site == "Seeley"],
                        mean_of_means$simulation[mean_of_means$site == "FivePoints"],
                        mean_of_means$simulation[mean_of_means$site == "McArthur"])
select_sites <- c("Seeley", "FivePoints", "McArthur")
lines(temperature ~ year, data = subset(averages, simulation == select_simulations[1]), lwd = 1, col = "orangered")
text(y = subset(averages, simulation == select_simulations[1] & year == 76)$temperature, x = 75 + 1, paste0(subset(averages, simulation == select_simulations[1] & year == 76)$site), pos = 4)
lines(temperature ~ year, data = subset(averages, simulation == select_simulations[2]), lwd = 1, col = "orange")
text(y = subset(averages, simulation == select_simulations[2] & year == 76)$temperature, x = 75 + 1, paste0(subset(averages, simulation == select_simulations[2] & year == 76)$site), pos = 4)
lines(temperature ~ year, data = subset(averages, simulation == select_simulations[3]), lwd = 1, col = "dodgerblue")
text(y = subset(averages, simulation == select_simulations[3] & year == 76)$temperature, x = 75 + 1, paste0(subset(averages, simulation == select_simulations[3] & year == 76)$site), pos = 4)
}; dev.off()

#### Plot weather station map ----

pdf("../figures/raw/fig_S10_weather_station_sites.pdf", height = 3, width = 3, pointsize = 3); {
plot(1, type = "n", xlim = c(-125, -115), ylim = c(10, 45), asp = 1.3, xlab = "", ylab = "", axes = FALSE)

map("worldHires", regions = c("USA", "Mexico"), xlim = c(-125, -110), ylim = c(10, 45),
    col = "gray90", fill = TRUE, lwd = 0.5, border = "gray40", add = TRUE)

map("state", regions = "california", add = TRUE, lwd = 0.5, col = "grey70", fill = TRUE, border = "black")

points(unique(dat$lon), unique(dat$lat), pch = 21, col = "gray40", lwd = 0.5, bg = "pink", cex = 1.5)
points(unique(dat$lon[dat$site == select_sites[1]]), unique(dat$lat[dat$site == select_sites[1]]), pch = 21, col = "gray40", lwd = 0.5, bg = "orangered", cex = 1.5)
points(unique(dat$lon[dat$site == select_sites[2]]), unique(dat$lat[dat$site == select_sites[2]]), pch = 21, col = "gray40", lwd = 0.5, bg = "orange", cex = 1.5)
points(unique(dat$lon[dat$site == select_sites[3]]), unique(dat$lat[dat$site == select_sites[3]]), pch = 21, col = "gray40", lwd = 0.5, bg = "dodgerblue", cex = 1.5)
}; dev.off()

#### Plot high-resolution temperature comparisons ----

subset_real_temperatures <- subset(dat, year(date.time) >= max(year(dat$date.time)) - 2 & month(date.time) == 6)
subset_simulated_temperatures <- subset(simulated_time_series, warming.rate == 0.04 & day >= 152 & day <= 152 + 29)
subset_simulated_temperatures <- subset_simulated_temperatures[seq(1, nrow(subset_simulated_temperatures), by = 2),]

## Real data
pdf("../figures/raw/fig_S10_real_temperature_data_high_res.pdf", height = 1.5, width = 3, pointsize = 3); {
plot(temp ~ as_datetime(date.time), subset(subset_real_temperatures, site == select_sites[1] & year(date.time) == max(year(subset_real_temperatures$date.time)) - 2), type = "l", ylab = "Temperature", xlab = "Date", ylim = c(0, 50), bty = "l", yaxt = "n", col = "orangered", lwd = 0.5)
lines(temp ~ as_datetime(date.time) %m-% years(1), subset(subset_real_temperatures, site == select_sites[2] & year(date.time) == max(year(subset_real_temperatures$date.time)) - 1), col = "orange", lwd = 0.5)
lines(temp ~ as_datetime(date.time) %m-% years(2), subset(subset_real_temperatures, site == select_sites[3] & year(date.time) == max(year(subset_real_temperatures$date.time))), col = "dodgerblue", lwd = 0.5)
axis(side = 4, at = c(5, 25, 45), las  = 2)
}; dev.off()

## Simulated data
pdf("../figures/raw/fig_S10_simulated_temperature_data_high_res.pdf", height = 1.5, width = 3, pointsize = 3); {
plot(temperature ~ time, subset(subset_simulated_temperatures, site == select_sites[1]), type = "l", ylab = "Temperature", ylim = c(0, 50), col = "orangered", xlab = "", xaxt = "n", yaxt = "n", lwd = 0.5)
lines(temperature ~ time, subset(subset_simulated_temperatures, site == select_sites[2]), type = "l", ylab = "Temperature", ylim = c(0, 50), col = "orange", xlab = "", xaxt = "n", yaxt = "n", lwd = 0.5)
lines(temperature ~ time, subset(subset_simulated_temperatures, site == select_sites[3]), type = "l", ylab = "Temperature", ylim = c(0, 50), col = "dodgerblue", xlab = "", xaxt = "n", yaxt = "n", lwd = 0.5)
axis(side = 4, at = c(5, 25, 45), las  = 2)
}; dev.off()

#### Validate simulations against observed data ----

## Simulate 5 years without warming. Cached to ../output/validation_time_series.txt so
## naive users never run this simulation (the file is provided on Zenodo and fetched by
## fetch_data.R). Preserve the 78-year series in main_time_series and restore it after the
## validation plots, so downstream scripts still receive the full series.
main_time_series <- simulated_time_series
if (!file.exists("../output/validation_time_series.txt")) {
  set.seed(599); simulated_time_series <- simulate_temperature(input_temperature = dat$temp, input_date_time = dat$date.time, input_site = dat$site, years = 5, n = 20, start_day = 182, end_day = 181)
  write.table(simulated_time_series, "../output/validation_time_series.txt", sep = "\t", row.names = F, quote = F)
} else {
  simulated_time_series <- read.delim("../output/validation_time_series.txt")
}

## Summary statistics
temperature_parameters <- as.data.frame(matrix(nrow = 20, ncol = 7))
colnames(temperature_parameters) <- c("site", "mean.true", "mean.sim", "var.true", "var.sim", "cor.acf", "cor.pacf")

for(i in 1:20){
  temperature_parameters$site[i] <- unique(dat$site)[i]
  temperature_parameters$mean.true[i] <- mean(subset(dat, site == unique(dat$site)[i])$temp)
  temperature_parameters$var.true[i] <- sd(subset(dat, site == unique(dat$site)[i])$temp)
  temperature_parameters$mean.sim[i] <- mean(subset(simulated_time_series, site == unique(dat$site)[i])$temperature)
  temperature_parameters$var.sim[i] <- sd(subset(simulated_time_series, site == unique(dat$site)[i])$temperature)
  temperature_parameters$cor.acf[i] <- cor(y = acf(aggregate(temp ~ as.Date(date.time), data = subset(dat, site == unique(dat$site)[i]), FUN = "mean"), plot = F, lag.max = 30)$acf,
                                           x = acf(aggregate(temperature ~ day, data = subset(simulated_time_series, site == unique(dat$site)[i]), FUN = "mean"), plot = F, lag.max = 30)$acf)
  temperature_parameters$cor.pacf[i] <- cor(y = pacf(subset(dat, site == unique(dat$site)[i])$temp, plot = F, lag.max = 24 * 3)$acf,
                                            x = pacf(subset(simulated_time_series, site == unique(dat$site)[i])$temperature, plot = F, lag.max = 24 * 3)$acf)
}

#### Validation plots ----

## Temperature distributions
pdf("../figures/raw/fig_S11_simu_vs_real_temp_dists.pdf", height = 1.5, width = 2, pointsize = 3); for(i in 1:20){
plot(density(subset(dat, site == unique(dat$site)[i])$temp), main = unique(dat$site)[i], ylab = "", yaxt = "n", xlab = "Temperature", xlim = c(-10, 55), lwd = 2)
lines(density(subset(simulated_time_series, site == unique(dat$site)[i])$temperature), lwd = 2, lty = 3, col = "red")
}; dev.off()

## Mean daily temperatures
pdf("../figures/raw/fig_S12_simu_vs_real_mean_daily_temps.pdf", height = 1.5, width = 2, pointsize = 3); {
daily_mean_obs <- aggregate(temp ~ day_of_year + site, data = dat[dat$day_of_year!=366,], FUN = mean)
daily_mean_sim <- aggregate(temperature ~ yday + site, data = simulated_time_series, FUN = mean)
for(i in 1:20){
  plot(temp ~ day_of_year, data = subset(daily_mean_obs, site == unique(dat$site)[i]), main = unique(dat$site)[i], xlab = "Day of year", ylab = "Temperature", ylim = c(-5, 40), xlim = c(1, 365), lwd = 1, type = "l")
    lines(temperature ~ yday, data = subset(daily_mean_sim, site == unique(dat$site)[i]), lwd = 1, lty = 3, col = "red")
}
}; dev.off()

## Mean hourly temperatures
pdf("../figures/raw/fig_S13_simu_vs_real_mean_hourly_temps.pdf", height = 1.5, width = 2, pointsize = 3); {
hourly_mean_obs <- aggregate(temp ~ hour + site, data = dat[dat$day_of_year!=366,], FUN = mean)
hourly_mean_sim <- aggregate(temperature ~ hour + site, data = simulated_time_series, FUN = mean)
for(i in 1:20){
    plot(temp ~ I(hour + 1), data = subset(hourly_mean_obs, site == unique(dat$site)[i]), main = unique(dat$site)[i], xlab = "Day of year", ylab = "Temperature", ylim = c(-5, 40), xlim = c(1, 24), lwd = 2, type = "l")
  lines(temperature ~ I(hour + 1), data = subset(hourly_mean_sim, site == unique(dat$site)[i]), lwd = 2, lty = 3, col = "red")
}
}; dev.off()

## Mean and variance comparison
pdf("../figures/raw/fig_S11_simu_vs_real_temp_parameters.pdf", height = 2, width = 2.5, pointsize = 3); {
  plot(temperature_parameters$mean.true ~ temperature_parameters$mean.sim, ylab = "Observed", xlab = "Simulated", cex = 2, pch = 21, bg = "orange", ylim = c(5, 20), xlim = c(5,20))
  points(temperature_parameters$var.true ~ temperature_parameters$var.sim, cex = 2, pch = 21, bg = "red3")
  abline(0, 1, lwd = 2, lty = 2)
}; dev.off()

## PACF correlation
pdf("../figures/raw/fig_S11_simu_vs_real_temp_pacf.pdf", height = 2, width = 1, pointsize = 3); {
  plot(temperature_parameters$cor.pacf ~ jitter(rep(1, 20), 1), ylab = "Observed", xlab = "Simulated", cex = 2, pch = 21, bg = "lightblue",xlim = c(0.5, 1.5), ylim = c(0, 1))
  abline(h = 1, lwd = 2, lty = 2)
  arrows(x0 = 0.95, x1 = 1.05, y0 = mean(temperature_parameters$cor.pacf), y1 = mean(temperature_parameters$cor.pacf), lwd = 2, code = 0)
}; dev.off()

## ACF correlation
pdf("../figures/raw/fig_S11_simu_vs_real_temp_acf.pdf", height = 2, width = 1, pointsize = 3); {
  plot(temperature_parameters$cor.acf ~ jitter(rep(1, 20), 1), ylab = "Observed", xlab = "Simulated", cex = 2, pch = 21, bg = "lightblue",xlim = c(0.5, 1.5), ylim = c(0, 1))
  abline(h = 1, lwd = 2, lty = 2)
  arrows(x0 = 0.95, x1 = 1.05, y0 = mean(temperature_parameters$cor.acf), y1 = mean(temperature_parameters$cor.acf), lwd = 2, code = 0)
}; dev.off()

## PACF time series
pdf("../figures/raw/fig_S11_simu_vs_real_temp_pacf_time_series.pdf", height = 1.5, width = 2, pointsize = 3); for(i in 1:20){
  plot(NA, ylim = c(-1, 1), xlim = c(0, 73), ylab = "PACF (hourly temperatures", xlab = "Lag (hours)", main = unique(dat$site)[i])
  abline(h = 0, lty = 2, lwd = 2)
  lines(pacf(subset(dat, site == unique(dat$site)[i])$temp, plot = F, lag.max = 24 * 3)$acf, lwd = 2)
  lines(pacf(subset(simulated_time_series, site == unique(dat$site)[i])$temperature, plot = F, lag.max = 24 * 3)$acf, col = "red", lwd = 2, lty = 3)
}; dev.off()

## ACF time series
pdf("../figures/raw/fig_S11_simu_vs_real_temp_acf_time_series.pdf", height = 1.5, width = 2, pointsize = 3); for(i in 1:20){
  plot(NA, ylim = c(-1, 1), xlim = c(0, 31), ylab = "ACF (daily average temperatures)", xlab = "Lag (days)", main = unique(dat$site)[i])
  abline(h = 0, lty = 2, lwd = 2)
  lines(acf(aggregate(temp ~ as.Date(date.time), data = subset(dat, site == unique(dat$site)[i]), FUN = "mean")$temp, plot = F, lag.max = 30)$acf, lwd = 2)
  lines(acf(aggregate(temperature ~ day, data = subset(simulated_time_series, site == unique(dat$site)[i]), FUN = "mean")$temperature, plot = F, lag.max = 30)$acf, col = "red", lwd = 2, lty = 3)
}; dev.off()

#### Restore the full 78-year series for downstream scripts ----
## The validation block above temporarily reassigned simulated_time_series; downstream
## scripts (results_20_sites*.R) need the 78-year series, so put it back.
simulated_time_series <- main_time_series
