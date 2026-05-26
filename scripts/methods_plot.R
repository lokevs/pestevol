#### Load packages, functions, and data ----
source("packages_functions_data.R")

#### Load data ----

TPC_parameters <- subset(rbind(read.table("../output/parameter_estimates_viab_fec.txt", header = T),
                               read.table("../output/parameter_estimates_dev_rate.txt", header = T),
                               read.table("../output/parameter_estimates_growth_rate.txt", header = T)),
                         effect == "main" & selection.regime == "anc" & origin == "bra" | effect == "main" & selection.regime == "anc" & trait == "viability")[,-c(3,6)]

plot_data <- subset(dat, selection.regime == "anc" & origin == "bra")

#### TPC fits to empirical data ----

## Fecundity
pdf("../figures/raw/methods_fec.pdf", height = 2, width = 2, pointsize = 3); {
set.seed(4596); plot(adult.offspring ~ jitter(temperature),
                     xlim = c(10, 50), ylim = c(0, 100), data = plot_data, pch = 21,
                     col = colorRampPalette(c("darkgrey", "white"))(6)[2],
                     bg = colorRampPalette(c("darkgrey", "white"))(6)[4],
                     ylab = "", xlab = "", lwd = 0.5)

lines(y = fec_function(temp = seq(10, 50, 0.1),
                       Topt = TPC_parameters$MAP[TPC_parameters$parameter == "Topt" & TPC_parameters$trait == "fecundity"],
                       fecmax = TPC_parameters$MAP[TPC_parameters$parameter == "fecmax" & TPC_parameters$trait == "fecundity"],
                       breadth = TPC_parameters$MAP[TPC_parameters$parameter == "breadth" & TPC_parameters$trait == "fecundity"],
                       fatness = TPC_parameters$MAP[TPC_parameters$parameter == "fatness" & TPC_parameters$trait == "fecundity"]),
      x = seq(10, 50, 0.1),
      col = "darkgrey", lwd = 2)

lines(y = 100 * invlogit(viab_function(temp = seq(10, 50, 0.1),
                                       a = TPC_parameters$MAP[TPC_parameters$parameter == "a" & TPC_parameters$trait == "viability"],
                                       h = TPC_parameters$MAP[TPC_parameters$parameter == "h" & TPC_parameters$trait == "viability"],
                                       k = TPC_parameters$MAP[TPC_parameters$parameter == "k" & TPC_parameters$trait == "viability"])),
      x = seq(10, 50, 0.1),
      col = colorRampPalette(c("darkgrey", 'black'))(6)[4], lwd = 2)

}; dev.off()

## Development rate
pdf("../figures/raw/methods_dev_rate.pdf", height = 2, width = 2, pointsize = 3); {
set.seed(4596); plot(dev.rate ~ jitter(temperature),
                     xlim = c(10, 50), ylim = c(0, 0.07), data = plot_data, pch = 21,
                     col = colorRampPalette(c("darkgrey", "white"))(6)[2],
                     bg = colorRampPalette(c("darkgrey", "white"))(6)[4],
                     ylab = "", xlab = "", lwd = 0.5)

lines(y = LRF(temp = seq(10, 50, 0.1),
              Tmin = TPC_parameters$MAP[TPC_parameters$parameter == "Tmin" & TPC_parameters$trait == "development_rate"],
              Tmax = TPC_parameters$MAP[TPC_parameters$parameter == "Tmax" & TPC_parameters$trait == "development_rate"],
              Topt = TPC_parameters$MAP[TPC_parameters$parameter == "Topt" & TPC_parameters$trait == "development_rate"],
              Ropt = TPC_parameters$MAP[TPC_parameters$parameter == "Ropt" & TPC_parameters$trait == "development_rate"]),
      x = seq(10, 50, 0.1),
      col = "darkgrey", lwd = 2)
}; dev.off()

## Growth rate
pdf("../figures/raw/methods_growth_rate.pdf", height = 2, width = 2, pointsize = 3); {
set.seed(4596); plot(growth.rate ~ jitter(temperature),
                     xlim = c(10, 50), ylim = c(0, 0.00011), data = plot_data[-75,], pch = 21,
                     col = colorRampPalette(c("darkgrey", "white"))(6)[2],
                     bg = colorRampPalette(c("darkgrey", "white"))(6)[4],
                     ylab = "", xlab = "", lwd = 0.5)

lines(y = LRF(temp = seq(10, 50, 0.1),
              Tmin = TPC_parameters$MAP[TPC_parameters$parameter == "Tmin" & TPC_parameters$trait == "growth_rate"],
              Tmax = TPC_parameters$MAP[TPC_parameters$parameter == "Tmax" & TPC_parameters$trait == "growth_rate"],
              Topt = TPC_parameters$MAP[TPC_parameters$parameter == "Topt" & TPC_parameters$trait == "growth_rate"],
              Ropt = TPC_parameters$MAP[TPC_parameters$parameter == "Ropt" & TPC_parameters$trait == "growth_rate"]),
      x = seq(10, 50, 0.1),
      col = "darkgrey", lwd = 2)
}; dev.off()

#### Hourly rates ----

## Oviposition rate
pdf("../figures/raw/methods_hourly_oviposition_rate.pdf", height = 1.5, width = 2, pointsize = 3); {
plot(y = oviposition_rate_function_anc_hourly(temp = seq(10, 50, 0.1), LRS = fec_function(temp = seq(10, 50, 0.1),
                                                                                                          Topt = TPC_parameters$MAP[TPC_parameters$parameter == "Topt" & TPC_parameters$trait == "fecundity"],
                                                                                                          fecmax = TPC_parameters$MAP[TPC_parameters$parameter == "fecmax" & TPC_parameters$trait == "fecundity"],
                                                                                                          breadth = TPC_parameters$MAP[TPC_parameters$parameter == "breadth" & TPC_parameters$trait == "fecundity"],
                                                                                                          fatness = TPC_parameters$MAP[TPC_parameters$parameter == "fatness" & TPC_parameters$trait == "fecundity"])),
                     x = seq(10, 50, 0.1),
                     xlim = c(10, 50), ylim = c(0, 1.5),
                     col = colorRampPalette(c("darkgrey", "white"))(6)[2],
                     bg = colorRampPalette(c("darkgrey", "white"))(6)[4],
                     ylab = "", xlab = "", lwd = 2, type = "l")
  }; dev.off()

## Viability rate
pdf("../figures/raw/methods_hourly_viability_rate.pdf", height = 1.5, width = 2, pointsize = 3); {
plot(y = viability_rate_hourly(temp = seq(10, 50, 0.1),
                                                                        Tmin = TPC_parameters$MAP[TPC_parameters$parameter == "Tmin" & TPC_parameters$trait == "development_rate"],
                                                                        Tmax = TPC_parameters$MAP[TPC_parameters$parameter == "Tmax" & TPC_parameters$trait == "development_rate"],
                                                                        Topt = TPC_parameters$MAP[TPC_parameters$parameter == "Topt" & TPC_parameters$trait == "development_rate"],
                                                                        Ropt = TPC_parameters$MAP[TPC_parameters$parameter == "Ropt" & TPC_parameters$trait == "development_rate"],
                                                                        a = TPC_parameters$MAP[TPC_parameters$parameter == "a" & TPC_parameters$trait == "viability"],
                                                                        h = TPC_parameters$MAP[TPC_parameters$parameter == "h" & TPC_parameters$trait == "viability"],
                                                                        k = TPC_parameters$MAP[TPC_parameters$parameter == "k" & TPC_parameters$trait == "viability"]),
                     x = seq(10, 50, 0.1),
                     xlim = c(10, 50), ylim = c(0.98, 1), ,
                     col = colorRampPalette(c("darkgrey", "white"))(6)[2],
                     bg = colorRampPalette(c("darkgrey", "white"))(6)[4],
                     ylab = "", xlab = "", lwd = 2, type = "l")
}; dev.off()

## Development rate (hourly)
pdf("../figures/raw/methods_hourly_dev_rate.pdf", height = 1.5, width = 2, pointsize = 3); {
plot(y = LRF_hourly(temp = seq(10, 50, 0.1),
                               Tmin = TPC_parameters$MAP[TPC_parameters$parameter == "Tmin" & TPC_parameters$trait == "development_rate"],
                               Tmax = TPC_parameters$MAP[TPC_parameters$parameter == "Tmax" & TPC_parameters$trait == "development_rate"],
                               Topt = TPC_parameters$MAP[TPC_parameters$parameter == "Topt" & TPC_parameters$trait == "development_rate"],
                               Ropt = TPC_parameters$MAP[TPC_parameters$parameter == "Ropt" & TPC_parameters$trait == "development_rate"]),
     x = seq(10, 50, 0.1),
     xlim = c(10, 50), ylim = c(0, 0.003), ,
     col = colorRampPalette(c("darkgrey", "white"))(6)[2],
     bg = colorRampPalette(c("darkgrey", "white"))(6)[4],
     ylab = "", xlab = "", lwd = 2, type = "l")
}; dev.off()

## Growth rate (hourly)
pdf("../figures/raw/methods_hourly_growth_rate.pdf", height = 1.5, width = 2, pointsize = 3); {
plot(y = LRF_hourly(temp = seq(10, 50, 0.1),
                    Tmin = TPC_parameters$MAP[TPC_parameters$parameter == "Tmin" & TPC_parameters$trait == "growth_rate"],
                    Tmax = TPC_parameters$MAP[TPC_parameters$parameter == "Tmax" & TPC_parameters$trait == "growth_rate"],
                    Topt = TPC_parameters$MAP[TPC_parameters$parameter == "Topt" & TPC_parameters$trait == "growth_rate"],
                    Ropt = TPC_parameters$MAP[TPC_parameters$parameter == "Ropt" & TPC_parameters$trait == "growth_rate"]),
     x = seq(10, 50, 0.1),
     xlim = c(10, 50), ylim = c(0, 0.000004), ,
     col = colorRampPalette(c("darkgrey", "white"))(6)[2],
     bg = colorRampPalette(c("darkgrey", "white"))(6)[4],
     ylab = "", xlab = "", lwd = 2, type = "l")
}; dev.off()

#### Yearly temperature trends ----

if (!exists("simulated_temperatures")) simulated_temperatures <- read.delim("../output/simulated_time_series.txt")
simulated_temperatures <- subset(simulated_temperatures, warming.rate == 0.04)

pdf("../figures/raw/methods_yearly_temps.pdf", height = 2, width = 3, pointsize = 3); {
plot(aggregate(temperature~year, data = subset(simulated_temperatures, site == "FivePoints"), "mean"), type = "l", bty = "l")
}; dev.off()

#### Simulation: year 2025 ----

temperatures_2025 <- subset(simulated_temperatures, year <= (min(simulated_temperatures$year) + 2))
temperatures_2025$hour.cumulative <- (temperatures_2025$day - 1) * 24 + temperatures_2025$hour
temperatures_2025$temp <- round(temperatures_2025$temp * 2)/2

# Expand missing origins with repeated values for each selection regime
TPC_parameters <- TPC_parameters %>%
  mutate(origin = ifelse(is.na(origin), list(c("bra", "ca", "yem")), origin)) %>%
  unnest(origin)
TPC_parameters <- as.data.frame(TPC_parameters)

# Extrinsic mortality parameters
lifespan <- 10
survival_rate <- exp(-(1 / (lifespan * 24)))

for(years in c(2)){
  for(wr in c(0.04)){

    for(si in "FivePoints"){

      temperatures <- subset(temperatures_2025, year <= years + 2 & year >= years -1 & site == si & warming.rate == wr)

      for(regime in c("anc")){
        for(orig in c("bra")){
          for(tr in c("all")){

            if(tr != "all" & regime != "hot"){
              next
            }

            regime_fec <- regime
            regime_viability <- regime
            regime_dev <- regime
            regime_growth <- regime

            if(tr == "fecundity" & regime == "hot"){
              regime_viability <- "anc"
              regime_dev <- "anc"
              regime_growth <- "anc"
            }
            if(tr == "viability" & regime == "hot"){
              regime_fec <- "anc"
              regime_dev <- "anc"
              regime_growth <- "anc"
            }
            if(tr == "development" & regime == "hot"){
              regime_fec <- "anc"
              regime_viability <- "anc"
              regime_growth <- "anc"
            }
            if(tr == "growth" & regime == "hot"){
              regime_fec <- "anc"
              regime_viability <- "anc"
              regime_dev <- "anc"
            }

            ## Initialize output
            simulated_fitness <- as.data.frame(matrix(nrow = 365, ncol = 15))
            colnames(simulated_fitness) <- c("year", "site", "warming.rate", "trait", "selection.regime", "origin", "day", "dev.time", "offspring.mass", "starting.mass", "viability", "offspring", "total.offspring.mass", "offspring.growth.rate", "fitness.rate")
            simulated_fitness$day <- 1:365
            simulated_fitness$year <- years+2023
            if(years == 2){
              simulated_fitness$warming.rate <- 0
            }
            else{
              simulated_fitness$warming.rate <- wr
            }
            simulated_fitness$trait <- tr
            simulated_fitness$selection.regime <- regime
            simulated_fitness$origin <- orig
            simulated_fitness$site <- si

            generate_TPCs(orig, regime_fec, regime_viability, regime_dev, regime_growth)

            ## Development time
            simulated_dev_time <- c()
            for(i in 1:(365 + 200)){
              temps <- temperatures[temperatures$day >= i + (years - 1) * 365,]$temp
              simulated_dev_time <- c(simulated_dev_time, which(cumsum(dev_rate(temps)) >= 1)[1] / 24)
            }
            simulated_fitness$dev.time <- simulated_dev_time[1:365]

            ## Growth
            simulated_mass <- c()
            for(i in 1:(365 + 200)){
              temps <- temperatures$temp[temperatures$day >= i + (years - 1) * 365][1:(simulated_dev_time[i] * 24)]
              growth <- sum(growth_rate(temps))
              simulated_mass <- c(simulated_mass, growth)
            }
            simulated_fitness$offspring.mass <- simulated_mass[1:365]

            ## Viability
            for(i in 1:365){
              temps <- temperatures$temp[temperatures$day >= i + (years - 1) * 365][1:(simulated_fitness$dev.time[i] * 24)]
              viability <- prod(viability_rate(temps))
              simulated_fitness$viability[i] <- viability
            }

            ## Starting mass
            for(i in 1:365){
              temps <- rev(temperatures$temp[temperatures$day < i + (years - 1) * 365])
              simulated_fitness$starting.mass[i] <-  sum(growth_rate(temps[1:which(cumsum(dev_rate(temps)) >= 1)[1]]))
            }

            ## Offspring production
            for(i in 1:365){
              temps <- temperatures$temp[temperatures$day >= i + (years - 1) * 365]
              fecmax_scaled <- fecmax * fecmax_scaling(simulated_fitness$starting.mass[i])

              max_counter <- min(which(cumsum(ovi_rate(temps)) >= fecmax_scaled)[1], which(survival_rate ^ (1:10000) <= 0.01)[1], na.rm = T)
              ovi_list <- ovi_rate(temps[1:max_counter]) * survival_rate ^ (1:max_counter)
              simulated_fitness$offspring[i] <- sum(ovi_list)
              daily_offspring <- colSums(matrix(c(ovi_list, rep(0, 24 - length(ovi_list) %% 24)), 24))

              # Apply viability to fitness calculations
              simulated_fitness$total.offspring.mass[i] <- sum(daily_offspring * simulated_mass[i:(i + length(daily_offspring) - 1)]) * simulated_fitness$viability[i]
              simulated_fitness$offspring.growth.rate[i] <- sum(daily_offspring * simulated_mass[i:(i + length(daily_offspring) - 1)] * (1 / (seq_along(daily_offspring) + simulated_dev_time[i:(i + length(daily_offspring) - 1)]))) * simulated_fitness$viability[i]
              simulated_fitness$fitness.rate[i] <- sum(daily_offspring * (1 / (seq_along(daily_offspring) + simulated_dev_time[i:(i + length(daily_offspring) - 1)]))) * simulated_fitness$viability[i]
            }
          }
        }
      }
    }
  }
}

fivepoints_bra_anc_2025 <- simulated_fitness

#### Simulation: year 2100 ----

temperatures_2100 <- subset(simulated_temperatures, year >= (max(simulated_temperatures$year) - 2))
temperatures_2100$hour.cumulative <- (temperatures_2100$day - 1) * 24 + temperatures_2100$hour
temperatures_2100$temp <- round(temperatures_2100$temp * 2)/2

for(years in c(77)){
  for(wr in c(0.04)){

    for(si in "FivePoints"){

      temperatures <- subset(temperatures_2100, year <= years + 2 & year >= years -1 & site == si & warming.rate == wr)

      for(regime in c("anc")){
        for(orig in c("bra")){
          for(tr in c("all")){

            if(tr != "all" & regime != "hot"){
              next
            }

            regime_fec <- regime
            regime_viability <- regime
            regime_dev <- regime
            regime_growth <- regime

            if(tr == "fecundity" & regime == "hot"){
              regime_viability <- "anc"
              regime_dev <- "anc"
              regime_growth <- "anc"
            }
            if(tr == "viability" & regime == "hot"){
              regime_fec <- "anc"
              regime_dev <- "anc"
              regime_growth <- "anc"
            }
            if(tr == "development" & regime == "hot"){
              regime_fec <- "anc"
              regime_viability <- "anc"
              regime_growth <- "anc"
            }
            if(tr == "growth" & regime == "hot"){
              regime_fec <- "anc"
              regime_viability <- "anc"
              regime_dev <- "anc"
            }

            ## Initialize output
            simulated_fitness <- as.data.frame(matrix(nrow = 365, ncol = 15))
            colnames(simulated_fitness) <- c("year", "site", "warming.rate", "trait", "selection.regime", "origin", "day", "dev.time", "offspring.mass", "starting.mass", "viability", "offspring", "total.offspring.mass", "offspring.growth.rate", "fitness.rate")
            simulated_fitness$day <- 1:365
            simulated_fitness$year <- years+2023
            if(years == 2){
              simulated_fitness$warming.rate <- 0
            }
            else{
              simulated_fitness$warming.rate <- wr
            }
            simulated_fitness$trait <- tr
            simulated_fitness$selection.regime <- regime
            simulated_fitness$origin <- orig
            simulated_fitness$site <- si

            generate_TPCs(orig, regime_fec, regime_viability, regime_dev, regime_growth)

            ## Development time
            simulated_dev_time <- c()
            for(i in 1:(365 + 200)){
              temps <- temperatures[temperatures$day >= i + (years - 1) * 365,]$temp
              simulated_dev_time <- c(simulated_dev_time, which(cumsum(dev_rate(temps)) >= 1)[1] / 24)
            }
            simulated_fitness$dev.time <- simulated_dev_time[1:365]

            ## Growth
            simulated_mass <- c()
            for(i in 1:(365 + 200)){
              temps <- temperatures$temp[temperatures$day >= i + (years - 1) * 365][1:(simulated_dev_time[i] * 24)]
              growth <- sum(growth_rate(temps))
              simulated_mass <- c(simulated_mass, growth)
            }
            simulated_fitness$offspring.mass <- simulated_mass[1:365]

            ## Viability
            for(i in 1:365){
              temps <- temperatures$temp[temperatures$day >= i + (years - 1) * 365][1:(simulated_fitness$dev.time[i] * 24)]
              viability <- prod(viability_rate(temps))
              simulated_fitness$viability[i] <- viability
            }

            ## Starting mass
            for(i in 1:365){
              temps <- rev(temperatures$temp[temperatures$day < i + (years - 1) * 365])
              simulated_fitness$starting.mass[i] <-  sum(growth_rate(temps[1:which(cumsum(dev_rate(temps)) >= 1)[1]]))
            }

            ## Offspring production
            for(i in 1:365){
              temps <- temperatures$temp[temperatures$day >= i + (years - 1) * 365]
              fecmax_scaled <- fecmax * fecmax_scaling(simulated_fitness$starting.mass[i])

              max_counter <- min(which(cumsum(ovi_rate(temps)) >= fecmax_scaled)[1], which(survival_rate ^ (1:10000) <= 0.01)[1], na.rm = T)
              ovi_list <- ovi_rate(temps[1:max_counter]) * survival_rate ^ (1:max_counter)
              simulated_fitness$offspring[i] <- sum(ovi_list)
              daily_offspring <- colSums(matrix(c(ovi_list, rep(0, 24 - length(ovi_list) %% 24)), 24))

              # Apply viability to fitness calculations
              simulated_fitness$total.offspring.mass[i] <- sum(daily_offspring * simulated_mass[i:(i + length(daily_offspring) - 1)]) * simulated_fitness$viability[i]
              simulated_fitness$offspring.growth.rate[i] <- sum(daily_offspring * simulated_mass[i:(i + length(daily_offspring) - 1)] * (1 / (seq_along(daily_offspring) + simulated_dev_time[i:(i + length(daily_offspring) - 1)]))) * simulated_fitness$viability[i]
              simulated_fitness$fitness.rate[i] <- sum(daily_offspring * (1 / (seq_along(daily_offspring) + simulated_dev_time[i:(i + length(daily_offspring) - 1)]))) * simulated_fitness$viability[i]
            }
          }
        }
      }
    }
  }
}

fivepoints_bra_anc_2100 <- simulated_fitness

#### Plot simulation results ----

## Trait predictions (2025 vs 2100)
pdf("../figures/raw/methods_trait_predictions.pdf", height = 1.5, width = 2, pointsize = 3)
plot(fivepoints_bra_anc_2025$dev.time ~ fivepoints_bra_anc_2025$day, ylab = "", main = "Development time (day)", ylim = c(0, 250),
     col = "darkgrey",xaxt = "n", xlab = "", lwd = 2, type = "l")
lines(fivepoints_bra_anc_2100$dev.time ~ fivepoints_bra_anc_2100$day, lty = 3,
     col = "darkgrey",xaxt = "n", xlab = "", lwd = 2, type = "l")

axis(1, at = c(1, 182, 365), labels = c("Jan 1", "Jul 1", "Dec 31"))

plot(fivepoints_bra_anc_2025$viability * 100 ~ fivepoints_bra_anc_2025$day, ylab = "", main = "Viability (%)", ylim = c(0, 100),
     col = "darkgrey",xaxt = "n", xlab = "", lwd = 2, type = "l")
lines(fivepoints_bra_anc_2100$viability * 100 ~ fivepoints_bra_anc_2100$day, lty = 3,
      col = "darkgrey",xaxt = "n", xlab = "", lwd = 2, type = "l")
axis(1, at = c(1, 182, 365), labels = c("Jan 1", "Jul 1", "Dec 31"))

plot(fivepoints_bra_anc_2025$offspring.mass ~ fivepoints_bra_anc_2025$day, ylab = "", main = "Offspring mass (mg)", ylim = c(0.0014, 0.00162), cex = 1,
     col = "darkgrey",xaxt = "n", xlab = "", lwd = 2, type = "l")
lines(fivepoints_bra_anc_2100$offspring.mass ~ fivepoints_bra_anc_2100$day, lty = 3,
      col = "darkgrey",xaxt = "n", xlab = "", lwd = 2, type = "l")
axis(1, at = c(1, 182, 365), labels = c("Jan 1", "Jul 1", "Dec 31"))

plot(fivepoints_bra_anc_2025$offspring ~ fivepoints_bra_anc_2025$day, ylab = "", main = "Adult offspring (high mort)", ylim = c(0, 100),
     col = "darkgrey",xaxt = "n", xlab = "", lwd = 2, type = "l")
lines(fivepoints_bra_anc_2100$offspring ~ fivepoints_bra_anc_2100$day, lty = 3,
      col = "darkgrey",xaxt = "n", xlab = "", lwd = 2, type = "l")
axis(1, at = c(1, 182, 365), labels = c("Jan 1", "Jul 1", "Dec 31"))
dev.off()

## Crop damage rate
pdf("../figures/raw/methods_crop_damage_rate_predictions.pdf", height = 2.5, width = 6, pointsize = 3)
plot(offspring.growth.rate ~ day, data = fivepoints_bra_anc_2025, ylab = "", main = "Crop damage rate", ylim = c(0, 0.004),
     col = "darkgrey",type = "l", xaxt = "n", lwd = 2)
lines(offspring.growth.rate ~ day, data = fivepoints_bra_anc_2100, lty = 3,
     col = "darkgrey",type = "l", lwd = 2)
axis(1, at = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335,365), labels = NA)
axis(1, at = c(15.5, 75, 136, 197, 259, 320), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"), tick = F)
dev.off()

## Fecmax-mass relationship
pdf("../figures/raw/methods_fecmax_mass_relationship.pdf", height = 1.2, width = 1.2, pointsize = 3)
fm <- lm(relative.adult.offspring ~ log(mean.weight) + selection.regime, data = subset(dat, temperature == 29 & selection.regime != "hotMcold"))
plot(relative.adult.offspring ~ (mean.weight), pch = 21, cex = 1.5, data = subset(dat, temperature == 29 & selection.regime == "anc" & origin == "bra"), bg = "darkgrey", bty = "l")
curve(coef(fm)[1] + coef(fm)[2] * log(x), 0.0005, 0.003, col = "darkgrey", lwd = 2, add = T)
dev.off()
