#### Setup ----
library(terra)
library(maps)
library(sf)
source("fitness_crop_damage_function_global.R") # for Fig 5B

## Load 12-band global rasters (1 deg, EPSG:4326)
## Bands per file: damage_{anc,hot,cold}, lambda_{...}, lambda_geo_{...}, repr_{...}
pre_data <- rast("../data/global_pres_mean.tif")
fut_data <- rast("../data/global_fut_mean.tif")

#### Shared layers ----

## Present (ancestral only)
dam_anc_pre <- pre_data[["damage_anc"]]
via_anc_pre <- pre_data[["lambda_geo_anc"]] >= 1

## Future ancestral baseline
dam_anc_fut <- fut_data[["damage_anc"]]
via_anc_fut <- fut_data[["lambda_geo_anc"]] >= 1

## Future per-pixel winning regime (anc > hot > cold ties)
l_anc  <- fut_data[["lambda_geo_anc"]]
l_hot  <- fut_data[["lambda_geo_hot"]]
l_cold <- fut_data[["lambda_geo_cold"]]
l_max  <- max(l_anc, l_hot, l_cold)
is_anc  <- l_anc == l_max
is_hot  <- (l_hot == l_max) & !is_anc
is_cold <- (l_cold == l_max) & !is_anc & !is_hot

via_best_fut <- l_max >= 1
dam_best_fut <- (fut_data[["damage_anc"]]  * is_anc) +
                (fut_data[["damage_hot"]]  * is_hot) +
                (fut_data[["damage_cold"]] * is_cold)

## Common map extent for panels
b_box <- ext(-180, 180, -50, 50)

#### Figure 6A top: present damage map ----

dam_pres_plot <- dam_anc_pre
dam_pres_plot[!via_anc_pre] <- NA
dam_pres_plot <- crop(dam_pres_plot, b_box)
cols_pres <- colorRampPalette(c("lightgreen", "purple3"))(100)
rng_pres <- minmax(dam_pres_plot)

pdf("../figures/raw/fig_6_present_damage_map.pdf", height = 2, width = 4, pointsize = 3)
plot(NA, xlim = c(-130, 180), ylim = c(-50, 50), axes = TRUE,
     main = "Current crop damage potential", xlab = "", ylab = "")
map("world", fill = TRUE, col = "grey93", border = "grey93", add = TRUE)
plot(dam_pres_plot, add = TRUE, col = cols_pres,
     legend = TRUE,
     plg = list(title = "Damage", title.cex = 0.8,
                at = c(rng_pres[1], rng_pres[2]), labels = c("Low", "High")))
abline(h = c(-40, -20, 0, 20, 40), lty = 2)
abline(v = c(-100, -50, 0, 50, 100, 150), lty = 2)

## Overlay all cowpea accessions (global, all source countries) as orange dots
accessions <- read.csv("../output/cowpea_accession_damage.csv")
points(accessions$lon, accessions$lat, pch = 16,
       col = "orange", cex = 0.2)

dev.off()

#### Figure 6A middle/bottom: fold-change maps ----

## log2 fold change in damage, NA where either present or future is non-viable
calc_fc <- function(d_future, d_present, v_future, v_present) {
  fc <- log2((d_future + 1e-9) / (d_present + 1e-9))
  fc[!(v_present & v_future)] <- NA
  fc
}

fc_warm  <- calc_fc(dam_anc_fut,  dam_anc_pre, via_anc_fut,  via_anc_pre)
fc_adapt <- calc_fc(dam_best_fut, dam_anc_pre, via_best_fut, via_anc_pre)

## Asymmetric blue-white-red palette centered on 0, shared across both maps
g_min <- min(minmax(fc_warm)[1], minmax(fc_adapt)[1], na.rm = TRUE)
g_max <- max(minmax(fc_warm)[2], minmax(fc_adapt)[2], na.rm = TRUE)
n_breaks <- 100
n_neg <- max(1, round(n_breaks * abs(g_min) / (g_max - g_min)))
n_pos <- max(1, n_breaks - n_neg)
cols_asym <- c(colorRampPalette(c("dodgerblue", "lightblue", "white"))(n_neg),
               colorRampPalette(c("white", "orange", "orangered"))(n_pos + 1)[-1])
at_vals <- seq(ceiling(g_min * 2)/2, floor(g_max * 2)/2, 0.5)

plot_fc_map <- function(fc_rast, v_num, v_den, main_title) {
  colonized <- !v_den & v_num # gained viability
  extinct   <- v_den & !v_num # lost viability
  colonized[colonized == 0] <- NA
  extinct[extinct == 0]     <- NA

  fc_rast   <- crop(fc_rast,   b_box)
  colonized <- crop(colonized, b_box)
  extinct   <- crop(extinct,   b_box)

  plot(NA, xlim = c(-130, 180), ylim = c(-50, 50), axes = TRUE,
       main = main_title, xlab = "", ylab = "")
  map("world", fill = TRUE, col = "grey93", border = "grey93", add = TRUE)
  plot(fc_rast, add = TRUE, col = cols_asym, range = c(g_min, g_max),
       legend = TRUE,
       plg = list(title = "log2(FC)", title.cex = 0.8,
                  at = at_vals, labels = at_vals))
  plot(colonized, add = TRUE, col = "darkred",  legend = FALSE)
  plot(extinct,   add = TRUE, col = "darkblue", legend = FALSE)
  abline(h = c(-40, -20, 0, 20, 40), lty = 2)
  abline(v = c(-100, -50, 0, 50, 100, 150), lty = 2)
  legend("bottomleft", legend = c("Colonization", "Extinction"),
         fill = c("darkred", "darkblue"), bty = "n", cex = 0.9)
}

pdf("../figures/raw/fig_6_fc_warming_map.pdf", height = 2, width = 4, pointsize = 3)
plot_fc_map(fc_warm, via_anc_fut, via_anc_pre, "Warming only")
dev.off()

pdf("../figures/raw/fig_6_fc_evolution_map.pdf", height = 2, width = 4, pointsize = 3)
plot_fc_map(fc_adapt, via_best_fut, via_anc_pre, "Warming and thermal adaptation")
dev.off()

#### Figure 6A right: latitudinal marginals ----

## Per-row mean across longitudes; NA only when the entire row is NA
lat_profile <- function(r) {
  lats <- yFromRow(r, 1:nrow(r))
  vals <- sapply(1:nrow(r), function(i) {
    v <- unlist(r[i, , drop = TRUE])
    if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  })
  data.frame(lat = lats, val = vals)
}

## Damage rasters with non-viable cells set to 0 (so they pull the mean down)
mb <- ext(-180, 180, -55, 55)
dam_pres_lat <- lat_profile(crop(dam_anc_pre  * via_anc_pre,  mb))
dam_warm_lat <- lat_profile(crop(dam_anc_fut  * via_anc_fut,  mb))
dam_best_lat <- lat_profile(crop(dam_best_fut * via_best_fut, mb))
occ_pres_lat <- lat_profile(crop(via_anc_pre,  mb))
occ_warm_lat <- lat_profile(crop(via_anc_fut,  mb))
occ_best_lat <- lat_profile(crop(via_best_fut, mb))

dam_xlim <- range(c(dam_pres_lat$val, dam_warm_lat$val, dam_best_lat$val), na.rm = TRUE)
ylims <- c(-50, 50)

pdf("../figures/raw/fig_6_marginals.pdf", height = 2, width = 1, pointsize = 3)

## Row 1: present
plot(dam_pres_lat$val, dam_pres_lat$lat, type = "l", col = "purple3", lwd = 1.5,
     ylim = ylims, xlim = dam_xlim, xaxt = "n", yaxt = "n", xlab = "Mean damage")
axis(1, at = dam_xlim, labels = signif(dam_xlim, 2))
abline(h = c(-40, -20, 0, 20, 40), lty = 2)
plot(occ_pres_lat$val, occ_pres_lat$lat, type = "l", col = "purple3", lwd = 1.5,
     ylim = ylims, xlim = c(0, 1), yaxt = "n", xlab = "P(λ ≥ 1)")
abline(h = c(-40, -20, 0, 20, 40), lty = 2)

## Row 2: warming only (present overlay for reference)
plot(dam_pres_lat$val, dam_pres_lat$lat, type = "l", col = "purple3", lwd = 1.5,
     ylim = ylims, xlim = dam_xlim, xaxt = "n", yaxt = "n", xlab = "Mean damage")
lines(dam_warm_lat$val, dam_warm_lat$lat, col = "orange", lwd = 1.5)
axis(1, at = dam_xlim, labels = signif(dam_xlim, 2))
abline(h = c(-40, -20, 0, 20, 40), lty = 2)
plot(occ_pres_lat$val, occ_pres_lat$lat, type = "l", col = "purple3", lwd = 1.5,
     ylim = ylims, xlim = c(0, 1), yaxt = "n", xlab = "P(λ ≥ 1)")
lines(occ_warm_lat$val, occ_warm_lat$lat, col = "orange", lwd = 1.5)
abline(h = c(-40, -20, 0, 20, 40), lty = 2)

## Row 3: warming + evolution
plot(dam_pres_lat$val, dam_pres_lat$lat, type = "l", col = "purple3", lwd = 1.5,
     ylim = ylims, xlim = dam_xlim, xaxt = "n", yaxt = "n", xlab = "Mean damage")
lines(dam_warm_lat$val, dam_warm_lat$lat, col = "orange", lwd = 1.5)
lines(dam_best_lat$val, dam_best_lat$lat, col = "red",    lwd = 1.5)
axis(1, at = dam_xlim, labels = signif(dam_xlim, 2))
abline(h = c(-40, -20, 0, 20, 40), lty = 2)
plot(occ_pres_lat$val, occ_pres_lat$lat, type = "l", col = "purple3", lwd = 1.5,
     ylim = ylims, xlim = c(0, 1), yaxt = "n", xlab = "P(λ ≥ 1)")
lines(occ_warm_lat$val, occ_warm_lat$lat, col = "orange", lwd = 1.5)
lines(occ_best_lat$val, occ_best_lat$lat, col = "red",    lwd = 1.5)
abline(h = c(-40, -20, 0, 20, 40), lty = 2)

dev.off()

#### Figure 5C: future winning regime (Robinson projection) ----

target_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

## Graticules at 90° lon, 30° lat
grat_proj <- project(vect(st_graticule(lon = seq(-180, 180, by = 90),
                                       lat = seq(-60, 60, by = 30))),
                     target_crs)

## Categorical raster: 1=anc, 2=hot, 3=cold; +3 if viable so 4-6 are the solid colors
win_cat  <- (is_anc * 1) + (is_hot * 2) + (is_cold * 3)
win_rast <- win_cat + (via_best_fut * 3)
win_rast[win_rast == 0] <- NA
win_rast <- project(crop(win_rast, ext(-180, 180, -60, 60)), target_crs, method = "near") # method="near" preserves integer codes

cols_solid <- c("darkgrey", "orangered", "dodgerblue") # anc / hot / cold
cols_trans <- adjustcolor(cols_solid, alpha.f = 0.5)

pdf("../figures/raw/fig_5_future_winners_map.pdf", height = 2.3, width = 3, pointsize = 3)
plot(NA, xlim = ext(win_rast)[1:2], ylim = ext(win_rast)[3:4],
     axes = FALSE, main = "Future winning thermal strategy", xlab = "", ylab = "")
plot(grat_proj, col = "grey80", lty = "dashed", add = TRUE)
plot(win_rast, add = TRUE, col = c(cols_trans, cols_solid),
     breaks = seq(0.5, 6.5, by = 1), legend = FALSE)
dev.off()

#### Figure 6B: environmental space, additive evolution benefit ----

mean_t <- rast("../data/global_tmean_2018_2024.tif")
sd_t   <- rast("../data/global_tsd_2018_2024.tif")

## Per-pixel future damage with non-viable regimes zeroed (extinct = no damage)
dam_anc_fut_v  <- as.vector(values(dam_anc_fut))
dam_best_fut_v <- as.vector(values(dam_best_fut))
via_anc_v      <- as.vector(values(via_anc_fut))
via_best_v     <- as.vector(values(via_best_fut))
dam_anc_fut_v[!via_anc_v]   <- 0
dam_best_fut_v[!via_best_v] <- 0

## Keep cells viable in at least one future scenario
keep <- via_anc_v | via_best_v
env_df <- na.omit(data.frame(
  tmean = as.vector(values(mean_t))[keep],
  tsd   = as.vector(values(sd_t))[keep],
  diff  = (dam_best_fut_v - dam_anc_fut_v)[keep]
))
env_df <- env_df[order(env_df$diff), ] # plot largest gains on top

## Symmetric palette so cells with no benefit sit at neutral grey
max_diff <- max(abs(env_df$diff), na.rm = TRUE)
breaks_sym <- seq(-max_diff, max_diff, length.out = 101)
pal_sym <- colorRampPalette(c("dodgerblue", "lightblue", "grey90", "orange", "orangered"))(100)
env_df$color <- pal_sym[as.numeric(cut(env_df$diff, breaks = breaks_sym, include.lowest = TRUE))]

pdf("../figures/raw/fig_6_env_space.pdf", height = 1.8, width = 1.6, pointsize = 3)
plot(env_df$tmean, env_df$tsd, xlim = c(9, 31),
     col = env_df$color, pch = 16, cex = 0.5,
     xlab = "Mean temperature (°C)", ylab = "Temperature SD",
     main = "Absolute impacts of evolution")
dev.off()

## Vertical gradient legend covering only the data range
rng <- range(env_df$diff, na.rm = TRUE)
idx <- as.numeric(cut(rng, breaks = breaks_sym, include.lowest = TRUE))
pal_used <- pal_sym[idx[1]:idx[2]]
n <- length(pal_used)
y_breaks <- seq(rng[1], rng[2], length.out = n + 1)

pdf("../figures/raw/fig_6_env_space_legend.pdf", height = 1.8, width = 0.7, pointsize = 3)
plot.new()
plot.window(xlim = c(0, 1), ylim = rng, xaxs = "i", yaxs = "i")
rect(0, y_breaks[-(n + 1)], 1, y_breaks[-1], col = pal_used, border = NA)
box(); axis(4, las = 1); mtext("Δ damage", side = 3, line = 0.3)
dev.off()

#### Figure 5B: constant-temperature winners ----

## At different constant temperatures, recompute lambda for each regime, plot bars
## colored by the winning regime at each step
temps <- seq(9.75, 45.25, 0.25)
tempdat <- data.frame(temp = temps, anc = NA, cold = NA, hot = NA)

for (i in seq_along(temps)) {
  tt <- calculate_fitness_crop_damage(rep(temps[i], 26280), plot = FALSE, lifespan_adult = 10)
  tempdat[i, "anc"]  <- tt[1, 4] # mean_lambda for anc
  tempdat[i, "cold"] <- tt[2, 4] # mean_lambda for cold
  tempdat[i, "hot"]  <- tt[3, 4] # mean_lambda for hot
}

lams <- as.matrix(tempdat[, c("anc", "cold", "hot")])
tempdat$max_lam  <- apply(lams, 1, max, na.rm = TRUE)
tempdat$winner   <- max.col(lams, ties.method = "first") # 1=anc, 2=cold, 3=hot
tempdat$bar_cols <- c("darkgrey", "dodgerblue", "red")[tempdat$winner]

pdf("../figures/raw/fig_5_constant_winners.pdf", height = 1.8, width = 1.8, pointsize = 3)
barplot(tempdat$max_lam, width = 1, space = 0,
        col = tempdat$bar_cols, border = NA, ylim = c(0, 1.2),
        yaxt = "n", xaxt = "n", xpd = FALSE,
        ylab = "Lambda", xlab = "Temperature")
axis(2, at = c(0, 1))
axis(1, at = (c(10, 20, 30, 40) - min(tempdat$temp)) / 0.25 + 0.5,
     labels = c(10, 20, 30, 40))
dev.off()
