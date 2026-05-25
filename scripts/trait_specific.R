#### Load packages ----
library(dplyr)

#### Data import ----

simu_dat <- subset(read.delim("../output/simu_dat_long.txt"), mortality == "med")

#### Crop damage: relative change ----

## Aggregate by trait, selection regime, origin, site, year, warming rate
simu_dat_aggregated_crop_damage <- aggregate(offspring.growth.rate ~ trait + selection.regime + origin + site + year + warming.rate, data = simu_dat, "mean")
simu_dat_aggregated_crop_damage$offspring.growth.change <- 0

## Calculate relative change from 2025 ancestor baseline
for(i in 1:nrow(simu_dat_aggregated_crop_damage)){
  subset_mask <- simu_dat_aggregated_crop_damage$selection.regime == "anc" &
    simu_dat_aggregated_crop_damage$site == simu_dat_aggregated_crop_damage$site[i] &
    simu_dat_aggregated_crop_damage$origin == simu_dat_aggregated_crop_damage$origin[i]

  anc_2025 <- simu_dat_aggregated_crop_damage$offspring.growth.rate[subset_mask & simu_dat_aggregated_crop_damage$warming.rate == 0]
  anc_current <- simu_dat_aggregated_crop_damage$offspring.growth.rate[subset_mask & simu_dat_aggregated_crop_damage$warming.rate == simu_dat_aggregated_crop_damage$warming.rate[i]]

  if(simu_dat_aggregated_crop_damage$selection.regime[i] == "anc"){
    simu_dat_aggregated_crop_damage$offspring.growth.change[i] <- round(((anc_current / anc_2025) - 1) * 100)
  } else {
    evo_current <- simu_dat_aggregated_crop_damage$offspring.growth.rate[i]
    simu_dat_aggregated_crop_damage$offspring.growth.change[i] <- round(((evo_current / anc_2025) - 1) * 100) - round(((anc_current / anc_2025) - 1) * 100)
  }
}

## Summarise crop damage across sites
simu_dat_aggregated_crop_damage_means <- simu_dat_aggregated_crop_damage %>%
  group_by(trait, selection.regime, origin, warming.rate) %>%
  summarise(
    min.value = min(offspring.growth.change, na.rm = TRUE),
    lower.quartile = quantile(offspring.growth.change, 0.25, na.rm = TRUE),
    mean.value = mean(offspring.growth.change, na.rm = TRUE),
    upper.quartile = quantile(offspring.growth.change, 0.75, na.rm = TRUE),
    max.value = max(offspring.growth.change, na.rm = TRUE),
    median.value = median(offspring.growth.change, na.rm = TRUE)
  ) %>% ungroup()

#### Fitness: relative change ----

## Aggregate by trait, selection regime, origin, site, year, warming rate
simu_dat_aggregated_fitness <- aggregate(fitness.rate ~ trait + selection.regime + origin + site + year + warming.rate, data = simu_dat, "mean")
simu_dat_aggregated_fitness$fitness.change <- 0

## Calculate relative change from 2025 ancestor baseline
for(i in 1:nrow(simu_dat_aggregated_fitness)){
  subset_mask <- simu_dat_aggregated_fitness$selection.regime == "anc" &
    simu_dat_aggregated_fitness$site == simu_dat_aggregated_fitness$site[i] &
    simu_dat_aggregated_fitness$origin == simu_dat_aggregated_fitness$origin[i]

  anc_2025 <- simu_dat_aggregated_fitness$fitness.rate[subset_mask & simu_dat_aggregated_fitness$warming.rate == 0]
  anc_current <- simu_dat_aggregated_fitness$fitness.rate[subset_mask & simu_dat_aggregated_fitness$warming.rate == simu_dat_aggregated_fitness$warming.rate[i]]

  if(simu_dat_aggregated_fitness$selection.regime[i] == "anc"){
    simu_dat_aggregated_fitness$fitness.change[i] <- round(((anc_current / anc_2025) - 1) * 100)
  } else {
    evo_current <- simu_dat_aggregated_fitness$fitness.rate[i]
    simu_dat_aggregated_fitness$fitness.change[i] <- round(((evo_current / anc_2025) - 1) * 100) - round(((anc_current / anc_2025) - 1) * 100)
  }
}

## Summarise fitness across sites
simu_dat_aggregated_fitness_means <- simu_dat_aggregated_fitness %>%
  group_by(trait, selection.regime, origin, warming.rate) %>%
  summarise(
    min.value = min(fitness.change, na.rm = TRUE),
    lower.quartile = quantile(fitness.change, 0.25, na.rm = TRUE),
    mean.value = mean(fitness.change, na.rm = TRUE),
    upper.quartile = quantile(fitness.change, 0.75, na.rm = TRUE),
    max.value = max(fitness.change, na.rm = TRUE),
    median.value = median(fitness.change, na.rm = TRUE)
  ) %>% ungroup()

#### Plotting ----

pdf("../figures/raw/fig_4_trait_specific_responses.pdf", height = 2.2, width = 1.25, pointsize = 3)

par(mfrow = c(1,1))
traits_to_plot <- c("all", "development", "fecundity", "growth", "viability")
origins <- c("bra", "yem", "ca")
lty_map <- c("bra" = 1, "yem" = 2, "ca" = 3)
pch_map <- c("bra" = 21, "yem" = 22, "ca" = 24)

add_custom_axes <- function() {
  axis(2, at = c(0, 100, 200))
  axis(1, at = c(0, 0.02, 0.04, 0.06))
}

## Crop damage: warming only (ancestral baseline)
plot(NA, xlim = c(-0.005, 0.065), ylim = c(-50, 200), main = "Warming only (ancestral baseline)",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
add_custom_axes()
abline(h = 0, lty = 2)

for(ori in origins){
  sub_dat <- subset(simu_dat_aggregated_crop_damage_means, origin == ori & selection.regime == "anc")
  lines(median.value ~ warming.rate, data = sub_dat, col = "darkgrey", lwd = 1.5, lty = lty_map[ori])
  points(median.value ~ warming.rate, data = sub_dat, col = colorRampPalette(c("darkgrey", "black"))(6)[4],
         bg = "darkgrey", lwd = 0.5, cex = 2, pch = pch_map[ori])
}

## Crop damage: evolution by trait
for(tr in traits_to_plot) {
  plot(NA, xlim = c(-0.005, 0.065), ylim = c(-50, 200), main = paste("Evolution of", tr),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  add_custom_axes()
  abline(h = 0, lty = 2)

  for(ori in origins){
    # Hot regime
    sub_hot <- subset(simu_dat_aggregated_crop_damage_means, trait == tr & origin == ori & selection.regime == "hot")
    lines(median.value ~ warming.rate, data = sub_hot, col = "orangered", lwd = 1.5, lty = lty_map[ori])
    points(median.value ~ warming.rate, data = sub_hot, col = colorRampPalette(c("orangered", "black"))(6)[4],
           bg = "orangered", lwd = 0.5, cex = 2, pch = pch_map[ori])

    # Cold regime
    sub_cold <- subset(simu_dat_aggregated_crop_damage_means, trait == tr & origin == ori & selection.regime == "cold")
    lines(median.value ~ warming.rate, data = sub_cold, col = "dodgerblue", lwd = 1.5, lty = lty_map[ori])
    points(median.value ~ warming.rate, data = sub_cold, col = colorRampPalette(c("dodgerblue", "black"))(6)[4],
           bg = "dodgerblue", lwd = 0.5, cex = 2, pch = pch_map[ori])
  }
}

## Fitness: warming only (ancestral baseline)
plot(NA, xlim = c(-0.005, 0.065), ylim = c(-50, 200), main = "Warming only (ancestral baseline)",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = 0, lty = 2)
add_custom_axes()

for(ori in origins){
  sub_dat <- subset(simu_dat_aggregated_fitness_means, origin == ori & selection.regime == "anc")
  lines(median.value ~ warming.rate, data = sub_dat, col = "darkgrey", lwd = 1.5, lty = lty_map[ori])
  points(median.value ~ warming.rate, data = sub_dat, col = colorRampPalette(c("darkgrey", "black"))(6)[4],
         bg = "darkgrey", lwd = 0.5, cex = 2, pch = pch_map[ori])
}

## Fitness: evolution by trait
for(tr in traits_to_plot) {
  plot(NA, xlim = c(-0.005, 0.065), ylim = c(-50, 200), main = paste("Evolution of", tr),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  add_custom_axes()
  abline(h = 0, lty = 2)

  for(ori in origins){
    # Hot regime
    sub_hot <- subset(simu_dat_aggregated_fitness_means, trait == tr & origin == ori & selection.regime == "hot")
    lines(median.value ~ warming.rate, data = sub_hot, col = "orangered", lwd = 1.5, lty = lty_map[ori])
    points(median.value ~ warming.rate, data = sub_hot, col = colorRampPalette(c("orangered", "black"))(6)[4],
           bg = "orangered", lwd = 0.5, cex = 2, pch = pch_map[ori])

    # Cold regime
    sub_cold <- subset(simu_dat_aggregated_fitness_means, trait == tr & origin == ori & selection.regime == "cold")
    lines(median.value ~ warming.rate, data = sub_cold, col = "dodgerblue", lwd = 1.5, lty = lty_map[ori])
    points(median.value ~ warming.rate, data = sub_cold, col = colorRampPalette(c("dodgerblue", "black"))(6)[4],
           bg = "dodgerblue", lwd = 0.5, cex = 2, pch = pch_map[ori])
  }
}

dev.off()
