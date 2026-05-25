source("packages_functions_data.R")

#### Read data ----

## TPC parameters
TPC_parameters <- subset(rbind(read.table("../output/parameter_estimates_viab_fec.txt", header = T),
                               read.table("../output/parameter_estimates_dev_rate.txt", header = T),
                               read.table("../output/parameter_estimates_growth_rate.txt", header = T)),
                         effect == "main" & selection.regime == "anc" & origin == "bra" | effect == "main" & selection.regime == "anc" & trait == "viability")[,-c(3,6)]

## Temperature data
temperatures <- read.delim("../data/climate_data_california.txt")
if (!exists("simulated_temperatures")) simulated_temperatures <- read.delim("../output/simulated_time_series.txt")
simulated_temperatures_04 <- subset(simulated_temperatures, warming.rate == 0.04)

#### Climate warming panel ----

pdf("../figures/raw/fig_1_intro_fig_climate_warming.pdf", height = 2, width = 2.5, pointsize = 3); {
  plot(aggregate(temperature~year, data = subset(simulated_temperatures_04, site == "FivePoints"), "mean"), type = "l", bty = "l")
  abline(lm(temperature~year, data = aggregate(temperature~year, data = subset(simulated_temperatures_04, site == "FivePoints"), "mean")), lty = 3)
}; dev.off()

#### Heat extremes panel ----

pdf("../figures/raw/fig_1_intro_fig_heat_extremes.pdf", height = 2, width = 2, pointsize = 3); {

## Temperature vectors
temp_cold <- temperatures$temp[temperatures$site == unique(temperatures$site)[3]][200:(200 + 24 * 7)]
temp_warm <- temp_cold + 3

## Time series
plot(temp_cold, type = "l", xlab = "Time", xaxt = "n", bty = "l",
     ylab = "Temperature", yaxt = "n", xlim = c(-10, 215), ylim = c(25, 48),
     col = "dodgerblue", lwd = 2)
lines(temp_warm, col = "orangered", lwd = 2)
axis(side = 2, at = c(30, 40), las = 2)

## Density estimates
dens_cold <- density(temp_cold, bw = 0.8)
dens_warm <- density(temp_warm, bw = 0.8)

scale_factor <- 250
x_offset <- 180

polygon(x = c(x_offset, x_offset + dens_cold$y * scale_factor, x_offset),
        y = c(dens_cold$x[1], dens_cold$x, dens_cold$x[length(dens_cold$x)]),
        col = rgb(30/255, 144/255, 255/255, 0.4),
        border = NA)

polygon(x = c(x_offset, x_offset + dens_warm$y * scale_factor, x_offset),
        y = c(dens_warm$x[1], dens_warm$x, dens_warm$x[length(dens_warm$x)]),
        col = rgb(255/255, 69/255, 0/255, 0.4),
        border = NA)

abline(h = 40, lty = 3)
abline(h = 35, lty = 3)

## Probability ratios of exceeding thermal thresholds
length(dens_warm$x[dens_warm$x >= 35]) / length(dens_cold$x[dens_cold$x >= 35])
length(dens_warm$x[dens_warm$x >= 40]) / length(dens_cold$x[dens_cold$x >= 40])

}; dev.off()

#### Performance extremes panel ----

pdf("../figures/raw/fig_1_intro_fig_performance_extremes.pdf", height = 2, width = 2, pointsize = 3); {

curve(LRF(x, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1), 14, 43.5, xlab = "Temperature", xaxt = "n", bty = "l",
     ylab = "Biological rate", yaxt = "n", xlim = c(16, 43), ylim = c(0, 1.1),
     col = "black", lwd = 2)
arrows(x0 = 30,
       x1 = 30,
       y0 = 0,
       y1 = LRF(30, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       lty = 1, lwd = 2, col = "orange", code = 0)
arrows(x0 = 0,
       x1 = 30,
       y0 = LRF(30, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       y1 = LRF(30, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       lty = 3, lwd = 2, col = "orange", code = 0)
arrows(x0 = 35,
       x1 = 35,
       y0 = 0,
       y1 = LRF(35, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       lty = 1, lwd = 2, col = "orangered", code = 0)
arrows(x0 = 0,
       x1 = 35,
       y0 = LRF(35, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       y1 = LRF(35, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       lty = 3, lwd = 2, col = "orangered", code = 0)
arrows(x0 = 40,
       x1 = 40,
       y0 = 0,
       y1 = LRF(40, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       lty = 1, lwd = 2, col = "darkred", code = 0)
arrows(x0 = 0,
       x1 = 40,
       y0 = LRF(40, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       y1 = LRF(40, Tmin = 17, Topt = 35, Tmax = 42, Ropt = 1),
       lty = 3, lwd = 2, col = "darkred", code = 0)

lines(temp_warm, col = "orangered", lwd = 2)
axis(side = 1, at = c(20, 30, 40), las = 1)

}; dev.off()

#### Daily means panel ----

pdf("../figures/raw/fig_1_intro_fig_daily_means.pdf", height = 2, width = 2, pointsize = 3); {

  temps <- c()
  for(i in 1:7){
  temps <- c(mean(temperatures$temp[temperatures$site == unique(temperatures$site)[3]][200:(199 + 24 * 7)][((i-1) * 24 + 1):((i) * 24)]), temps)
  }

  plot(y = rep(temps, each = 24), x = 1:(24*7), type = "l", xlab = "Time", xaxt = "n", bty = "l",
       ylab = "Temperature", yaxt = "n", xlim = c(-10, 215), ylim = c(25, 48),
       col = "dodgerblue", lwd = 2)
  axis(side = 2, at = c(30, 40), las = 2)

}; dev.off()
