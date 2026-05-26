source("packages_functions_data.R")

#### Plot raw data ----

pdf("../figures/raw/fig_2_raw_data.pdf", height = 2, width = 1.5, pointsize = 3); {
  set.seed(5285)

  ## Development rate
  for(g in c("bra", "ca", "yem")) {

    plot(NA,
         ylab = "Development rate", xlab = "Temperature", xaxt = "n", yaxt = "n",
         ylim = c(0, 0.065), xlim = c(14, 40),
         main = g)
    abline(h = c(0, 0.025, 0.05), lty = 3, col = "grey")

    for(i in c("anc", "cold", "hot")){
      tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

      points(jitter(dev.rate) ~ jitter(temperature),
             data = tempdat,
             pch = unique(tempdat$point), cex = 1,
             lwd = 0.5,
             bg = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[4],
             col = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[2])
    }

    groupmeans_regime <- aggregate(dev.rate ~ replicate + temperature + origin + colour + point, data = subset(dat, origin == g & both.sexes == 1 & adult.offspring != 0), "geomean")

    # Lines connecting group means
    for(re in unique(groupmeans_regime$replicate)){
      tempdat <- subset(groupmeans_regime, replicate == re)

      lines(dev.rate ~ temperature,
            data = tempdat,
            lwd = 1,
            lty = unique(tempdat$linetype),
            col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
    }

    # Group mean points
    for(re in unique(groupmeans_regime$replicate)){
      tempdat <- subset(groupmeans_regime, replicate == re)

      points(dev.rate ~ temperature,
             data = tempdat,
             pch = unique(tempdat$point), cex = 2,
             lwd = 0.5,
             lty = 3,
             bg = unique(tempdat$colour),
             col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
    }

    axis(2, at = c(0, 0.025, 0.05), labels = c(0, "", 0.05))
    axis(1, at = c(15, 25, 35))
  }

  legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n")

  ## Adult offspring
  for(g in c("bra", "ca", "yem")) {

    plot(NA,
         ylab = "Adult offspring", xlab = "Temperature", xaxt = "n", yaxt = "n",
         ylim = c(0, 120), xlim = c(14, 40),
         main = g)
    abline(h = c(0, 30, 60, 90, 120), lty = 3, col = "grey")

    for(i in c("anc", "cold", "hot")){
      tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

      points(jitter(adult.offspring) ~ jitter(temperature),
             data = tempdat,
             pch = unique(tempdat$point), cex = 1,
             lwd = 0.5,
             bg = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[4],
             col = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[2])
    }

    groupmeans_regime <- aggregate(adult.offspring ~ replicate + temperature + origin + colour + point, data = subset(dat, origin == g & both.sexes == 1 & adult.offspring != 0), "geomean")

    # Lines connecting group means
    for(re in unique(groupmeans_regime$replicate)){
      tempdat <- subset(groupmeans_regime, replicate == re)

      lines(adult.offspring ~ temperature,
            data = tempdat,
            lwd = 1,
            lty = unique(tempdat$linetype),
            col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
    }

    # Group mean points
    for(re in unique(groupmeans_regime$replicate)){
      tempdat <- subset(groupmeans_regime, replicate == re)

      points(adult.offspring ~ temperature,
             data = tempdat,
             pch = unique(tempdat$point), cex = 2,
             lwd = 0.5,
             lty = 3,
             bg = unique(tempdat$colour),
             col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
    }

    axis(2, at = c(0, 60, 120))
    axis(1, at = c(15, 25, 35))
  }

  ## Growth rate
  for(g in c("bra", "ca", "yem")) {

    plot(NA,
         ylab = "Growth rate", xlab = "Temperature", xaxt = "n", yaxt = "n",
         ylim = c(0, 0.15/1000), xlim = c(14, 40),
         main = g)
    abline(h = c(0, 0.05, 0.1, 0.15)/1000, lty = 3, col = "grey")

    for(i in c("anc", "cold", "hot")){
      tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

      points(jitter(growth.rate) ~ jitter(temperature),
             data = tempdat,
             pch = unique(tempdat$point), cex = 1,
             lwd = 0.5,
             bg = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[4],
             col = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[2])
    }

    groupmeans_regime <- aggregate(growth.rate ~ replicate + temperature + origin + colour + point, data = subset(dat, origin == g & both.sexes == 1 & adult.offspring != 0), "geomean")

    # Lines connecting group means
    for(re in unique(groupmeans_regime$replicate)){
      tempdat <- subset(groupmeans_regime, replicate == re)

      lines(growth.rate ~ temperature,
            data = tempdat,
            lwd = 1,
            lty = unique(tempdat$linetype),
            col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
    }

    # Group mean points
    for(re in unique(groupmeans_regime$replicate)){
      tempdat <- subset(groupmeans_regime, replicate == re)

      points(growth.rate ~ temperature,
             data = tempdat,
             pch = unique(tempdat$point), cex = 2,
             lwd = 0.5,
             lty = 3,
             bg = unique(tempdat$colour),
             col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
    }

    axis(2, at = c(0, 0.05, 0.1, 0.15)/1000, labels = c("0", "", "0.1", ""))
    axis(1, at = c(15, 25, 35))
  }
}; dev.off()

#### Plot relative data ----

pdf("../figures/raw/fig_2_raw_relative_data.pdf", height = 1.2, width = 1.5, pointsize = 3); {

  ## Relative development rate
  tempdat <- aggregate(dev.rate ~ selection.regime + origin + temperature + replicate + colour + point, data = subset(dat, both.sexes == 1),  "mean")
  tempdat$dev.rate <- as.numeric(tempdat$dev.rate)
  tempdat$temperature <- as.numeric(tempdat$temperature)
  tempdat$point <- as.integer(tempdat$point)
  tempdat$dev.rate.relative <- 0

  for(i in 1:nrow(tempdat)){
    tempdat$dev.rate.relative[i] <- tempdat$dev.rate[i] / max(tempdat$dev.rate[tempdat$origin == tempdat$origin[i] & tempdat$temperature == tempdat$temperature[i]])
  }

  for(g in c("bra", "ca", "yem")) {
    plot(NA,
         ylab = "", xlab = "", xaxt = "n", yaxt = "n",
         ylim = c(-0.05 + min(tempdat$dev.rate.relative, na.rm = T), 1.05),
         xlim = c(14, 40),
         main = g)
    axis(2,  round(c(min(tempdat$dev.rate.relative, na.rm = T), 1), 2))

    for(re in unique(tempdat$replicate)){
      tempdat2 <- subset(tempdat, replicate == re & origin == g)
      tempdat2 <- tempdat2[order(tempdat2$temperature),]

      lines(dev.rate.relative ~ temperature,
            data = tempdat2,
            lwd = 2,
            col = alpha(unique(tempdat2$colour), 0.5))
    }
  }

  ## Relative adult offspring
  tempdat <- aggregate(adult.offspring ~ selection.regime + origin + temperature + replicate + colour + point, data = subset(dat, both.sexes == 1),  "mean")
  tempdat$adult.offspring.relative <- 0

  for(i in 1:nrow(tempdat)){
    tempdat$adult.offspring.relative[i] <- tempdat$adult.offspring[i] / max(tempdat$adult.offspring[tempdat$origin == tempdat$origin[i] & tempdat$temperature == tempdat$temperature[i]])
  }

  for(g in c("bra", "ca", "yem")) {
    plot(NA,
         ylab = "", xlab = "", xaxt = "n", yaxt = "n",
         ylim = c(-0.05 + min(tempdat$adult.offspring.relative, na.rm = T), 1.05),
         xlim = c(14, 40),
         main = g)
    axis(2,  round(c(min(tempdat$adult.offspring.relative, na.rm = T), 1), 2))

    for(re in unique(tempdat$replicate)){
      tempdat2 <- subset(tempdat, replicate == re & origin == g)
      tempdat2 <- tempdat2[order(tempdat2$temperature),]

      lines(adult.offspring.relative ~ temperature,
            data = tempdat2,
            lwd = 2,
            col = alpha(unique(tempdat2$colour), 0.5))
    }
  }

  ## Relative growth rate
  tempdat <- aggregate(growth.rate ~ selection.regime + origin + temperature + replicate + colour + point, data = subset(dat, both.sexes == 1),  "mean")
  tempdat$growth.rate <- as.numeric(tempdat$growth.rate)
  tempdat$temperature <- as.numeric(tempdat$temperature)
  tempdat$point <- as.integer(tempdat$point)
  tempdat$growth.rate.relative <- 0

  for(i in 1:nrow(tempdat)){
    tempdat$growth.rate.relative[i] <- tempdat$growth.rate[i] / max(tempdat$growth.rate[tempdat$origin == tempdat$origin[i] & tempdat$temperature == tempdat$temperature[i]])
  }

  for(g in c("bra", "ca", "yem")) {
    plot(NA,
         ylab = "", xlab = "", xaxt = "n", yaxt = "n",
         ylim = c(-0.05 + min(tempdat$growth.rate.relative, na.rm = T), 1.05),
         xlim = c(14, 40),
         main = g)
    axis(2,  round(c(min(tempdat$growth.rate.relative, na.rm = T), 1), 2))

    for(re in unique(tempdat$replicate)){
      tempdat2 <- subset(tempdat, replicate == re & origin == g)
      tempdat2 <- tempdat2[order(tempdat2$temperature),]

      lines(growth.rate.relative ~ temperature,
            data = tempdat2,
            lwd = 2,
            col = alpha(unique(tempdat2$colour), 0.5))
    }
  }
}; dev.off()

#### Origin map ----

world_data <- map_data("world")

# Filter out small island polygons
world_data_grouped <- world_data %>%
  group_by(group) %>%
  mutate(n_points = n()) %>%
  ungroup()

min_points <- 100
world_data_filtered <- world_data_grouped %>%
  filter(n_points >= min_points)

highlight_countries <- c("USA", "Brazil", "Yemen")

highlight_data <- world_data_filtered %>%
  filter(region %in% highlight_countries)

gg_map_custom <- ggplot() +
  geom_polygon(
    data = world_data_filtered,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey92",
    linewidth = 0.2
  ) +
  geom_polygon(
    data = highlight_data,
    aes(x = long, y = lat, group = group),
    fill = "lightblue",
    color = NA
  ) +
  coord_quickmap(xlim = c(-130, 60), ylim = c(-40, 60)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

gg_map_custom

ggsave(
  filename = "../figures/raw/fig_2_map.pdf",
  plot = gg_map_custom,
  width = 1.5,
  height = 1.2,
  units = "in"
)
