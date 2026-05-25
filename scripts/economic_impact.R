#### Setup ----

library(terra)
library(dplyr)
library(genesysr)
library(sf)
library(rnaturalearth)
library(ggplot2)

#### Cowpea production countries ----

most_farmed <- read.delim("../data/cowpea_production.txt")

iso3_lookup <- c(
  "Nigeria" = "NGA", "Brazil" = "BRA", "Niger" = "NER",
  "Burkina Faso" = "BFA", "Ghana" = "GHA",
  "Cameroon" = "CMR", "Kenya" = "KEN", "Mali" = "MLI",
  "Mozambique" = "MOZ", "Myanmar" = "MMR", "Sudan" = "SDN",
  "Tanzania" = "TZA",
  "Bosnia and Herzegovina" = "BIH", "Botswana" = "BWA",
  "Chad" = "TCD", "China" = "CHN", "Croatia" = "HRV",
  "DRC" = "COD", "Egypt" = "EGY", "Guinea" = "GIN",
  "Guyana" = "GUY", "Hungary" = "HUN", "Iraq" = "IRQ",
  "Japan" = "JPN", "Madagascar" = "MDG", "Malawi" = "MWI",
  "Mauritania" = "MRT", "Namibia" = "NAM", "North Macedonia" = "MKD",
  "Peru" = "PER", "Philippines" = "PHL", "Senegal" = "SEN",
  "Serbia" = "SRB", "Sierra Leone" = "SLE", "South Africa" = "ZAF",
  "South Sudan" = "SSD", "Sri Lanka" = "LKA", "Togo" = "TGO",
  "Uganda" = "UGA", "USA" = "USA",
  "India" = "IND" # not in Boukar et al.; no quantified crop production
)
country_codes <- unname(iso3_lookup[most_farmed$Country])
stopifnot(!any(is.na(country_codes)))

#### Load rasters and apply extinction mask ----

## 12-band global rasters; bands include damage_{anc,hot,cold} and lambda_geo_{...}
pres <- rast("../data/global_pres_mean.tif")
fut  <- rast("../data/global_fut_mean.tif")

## Per-pixel damage with non-persistent regimes set to 0 (lambda_geo < 1 = extinction)
mask_extinct <- function(damage, lambda_geo) ifel(lambda_geo >= 1, damage, 0)

damage_pres_m    <- mask_extinct(pres$damage_anc, pres$lambda_geo_anc)
damage_fut_anc_m <- mask_extinct(fut$damage_anc,  fut$lambda_geo_anc)

## Future-with-evolution: take the genotype with highest lambda_geo per pixel;
## if even the best is < 1, the cell is extinct (damage = 0)
lam_stack  <- c(fut$lambda_geo_anc, fut$lambda_geo_hot, fut$lambda_geo_cold)
dmg_stack  <- c(fut$damage_anc,     fut$damage_hot,     fut$damage_cold)
winner     <- which.max(lam_stack)
dmg_winner <- ifel(winner == 1, dmg_stack[[1]],
                   ifel(winner == 2, dmg_stack[[2]], dmg_stack[[3]]))
damage_fut_evo_m <- ifel(max(lam_stack) >= 1, dmg_winner, 0)

#### Cowpea accessions from Genesys (cached) ----

## Server-side countryOfOrigin filter is unreliable across API versions;
## download all Vigna and filter by ORIGCTY locally. Result is cached so
## repeat runs do not need OAuth.

accession_cache <- "../output/cowpea_accessions_geo.csv"

if (file.exists(accession_cache)) {
  cowpea_geo <- read.csv(accession_cache)
} else {
  user_login() # opens browser for OAuth; once per session

  cowpea <- get_accessions(
    filters = list(taxonomy = list(genus = list("Vigna"))),
    fields  = list("INSTCODE", "ACCENUMB", "DOI",
                   "GENUS", "SPECIES", "SUBTAXA",
                   "ORIGCTY",
                   "DECLATITUDE", "DECLONGITUDE", "ELEVATION",
                   "COLLSITE", "COLLDATE", "SAMPSTAT", "COLLSRC")
  )

  cowpea_geo <- cowpea %>%
    rename(
      acceNumb = ACCENUMB, instCode = INSTCODE, doi = DOI,
      genus    = GENUS,    species  = SPECIES,
      origcty  = ORIGCTY,
      lat      = DECLATITUDE, lon = DECLONGITUDE, elev = ELEVATION
    ) %>%
    filter(species == "unguiculata",
           !is.na(lat), !is.na(lon),
           COLLSRC %in% 20:28) # farm or cultivated habitat (MCPD); kept globally

  # Country filtering is deferred: Fig 6A uses the global set; Fig 6C aggregates
  # only the producer countries + India.
  write.csv(cowpea_geo, accession_cache, row.names = FALSE)
}

## Collapse accessions sharing exact coordinates (treat as one sampling location)
cowpea_geo <- cowpea_geo %>% distinct(lat, lon, .keep_all = TRUE)

#### Per-accession damage ----

## Damage values are extracted from the extinction-masked rasters so accessions
## in non-viable pixels contribute 0 (matches the convention used everywhere else)
pts <- vect(cowpea_geo, geom = c("lon", "lat"), crs = "EPSG:4326")
cowpea_geo$damage_present       <- terra::extract(damage_pres_m,    pts, ID = FALSE)[, 1]
cowpea_geo$damage_fut_no_evol   <- terra::extract(damage_fut_anc_m, pts, ID = FALSE)[, 1]
cowpea_geo$damage_fut_with_evol <- terra::extract(damage_fut_evo_m, pts, ID = FALSE)[, 1]

## Winning genotype is identified from the raw (unmasked) lambda values per regime
fut_vals   <- terra::extract(fut, pts, ID = FALSE)
lambda_mat <- as.matrix(fut_vals[, c("lambda_geo_anc", "lambda_geo_hot", "lambda_geo_cold")])
w_acc      <- max.col(lambda_mat, ties.method = "first") # 1=anc, 2=hot, 3=cold
cowpea_geo$winning_genotype <- c("anc", "hot", "cold")[w_acc]

## Accessions outside the raster domain
na_rows <- rowSums(is.na(lambda_mat)) == 3
cowpea_geo$winning_genotype[na_rows]     <- NA
cowpea_geo$damage_fut_with_evol[na_rows] <- NA

cowpea_geo$delta_no_evol   <- cowpea_geo$damage_fut_no_evol   - cowpea_geo$damage_present
cowpea_geo$delta_with_evol <- cowpea_geo$damage_fut_with_evol - cowpea_geo$damage_present

write.csv(cowpea_geo, "../output/cowpea_accession_damage.csv", row.names = FALSE)

#### Country-level accession-weighted means ----

## Each accession is one observation at its pixel's damage value; averaging across
## accessions per country naturally weights each pixel by accession density.
## Countries with only one sampling location have zero within-country variation,
## so they are dropped; the minimum kept is n = 2 (also excludes MMR, MOZ).
country_damage <- cowpea_geo %>%
  filter(!is.na(damage_present),
         !is.na(damage_fut_no_evol),
         !is.na(damage_fut_with_evol)) %>%
  group_by(origcty) %>%
  summarise(
    n_accessions         = n(),
    damage_present       = mean(damage_present),
    damage_fut_no_evol   = mean(damage_fut_no_evol),
    damage_fut_with_evol = mean(damage_fut_with_evol),
    .groups = "drop"
  ) %>%
  filter(n_accessions >= 2) %>%
  rename(iso3 = origcty) %>%
  mutate(country         = names(iso3_lookup)[match(iso3, iso3_lookup)],
         delta_no_evol   = damage_fut_no_evol   - damage_present,
         delta_with_evol = damage_fut_with_evol - damage_present) %>%
  select(iso3, country, n_accessions,
         damage_present, damage_fut_no_evol, damage_fut_with_evol,
         delta_no_evol,  delta_with_evol) %>%
  arrange(desc(damage_fut_no_evol))

write.csv(country_damage, "../output/cowpea_country_damage.csv", row.names = FALSE)

#### Figure 6A overlay: cowpea accessions ----

world            <- ne_countries(scale = "medium", returnclass = "sf")
cowpea_countries <- world[world$iso_a3 %in% country_codes, ]

map_accessions <- ggplot() +
  geom_sf(data = world,            fill = NA,       colour = "grey90", linewidth = 0.1) +
  geom_sf(data = cowpea_countries, fill = "grey96", colour = "grey60", linewidth = 0.3) +
  geom_point(data = cowpea_geo,
             aes(x = lon, y = lat, colour = winning_genotype),
             size = 0.7, alpha = 0.6) +
  scale_colour_manual(
    values   = c(anc = "#7F7F7F", hot = "#D55E00", cold = "#0072B2"),
    na.value = "black",
    name     = "Winning genotype\n(future climate)",
    labels   = c(anc = "Ancestral", hot = "Hot-adapted", cold = "Cold-adapted")
  ) +
  coord_sf(xlim = c(-130, 180), ylim = c(-50, 50), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom")

ggsave("../figures/raw/fig_6_cowpea_accessions.pdf",
       plot = map_accessions, width = 4, height = 2)

#### Figure 6C: country-level damage barplot ----

## Order countries by production tier, then by future-with-evol damage within each tier.
## Tier join uses iso3 (not country name). India has no production tier and is
## appended at the right with a pink bar as a non-quantified reference.
tier_levels <- unique(most_farmed[[2]])
country_damage$Tier <- factor(
  most_farmed[[2]][match(country_damage$iso3, country_codes)],
  levels = tier_levels
)
dat <- country_damage[!is.na(country_damage$Tier), ]
dat <- dat[order(dat$Tier, -dat$damage_fut_with_evol), ]

## Append India at the rightmost position, if present after filtering
india_row <- country_damage[country_damage$iso3 == "IND", ]
dat <- rbind(dat, india_row)

## Bar color encodes production tier: darkest = top producer, lightest = smallest
tier_colors <- colorRampPalette(c("grey30", "grey85"))(length(tier_levels))
names(tier_colors) <- tier_levels
bar_colors <- tier_colors[as.character(dat$Tier)]
bar_colors[dat$iso3 == "IND"] <- "#FF99CC" # India: non-quantified reference

## Slightly larger gap before India to visually separate it from producers
spaces <- rep(0.2, nrow(dat))
if (any(dat$iso3 == "IND")) spaces[which(dat$iso3 == "IND")] <- 0.8

pdf("../figures/raw/fig_6_country_bars.pdf", height = 2, width = 3, pointsize = 3)
bp <- barplot(dat$damage_present, names.arg = dat$country, las = 2,
              space = spaces, col = bar_colors, border = NA, ylab = "Damage Potential",
              ylim = c(0, max(dat[, c("damage_present", "damage_fut_no_evol",
                                      "damage_fut_with_evol")], na.rm = TRUE) * 1.05),
              cex.names = 0.8)
arrows(x0 = bp, y0 = dat$damage_fut_no_evol,
       x1 = bp, y1 = dat$damage_fut_with_evol,
       col = "#0072B2", length = 0.05, lwd = 2)
arrows(x0 = bp, y0 = dat$damage_present,
       x1 = bp, y1 = dat$damage_fut_no_evol,
       col = "#D55E00", length = 0.05, lwd = 2)
legend("topright",
       legend = c("Present", "Trajectory: Warming", "Trajectory: Evolution"),
       col    = c("grey85",  "#D55E00",             "#0072B2"),
       lwd    = c(8, 2, 2), bty = "n")
dev.off()

#### Figure 6C inset: present-damage correlations ----

dat$only_evol <- dat$damage_fut_with_evol - dat$damage_fut_no_evol
dat$only_warm <- dat$damage_fut_no_evol   - dat$damage_present

pdf("../figures/raw/fig_6_country_correlations.pdf", height = 1.2, width = 1, pointsize = 3)

plot(only_warm ~ damage_present, dat,
     pch = 21, cex = 1, bg = "grey", lwd = 0.25,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     xlim = c(min(damage_present)-0.0001, max(damage_present)+0.0001),
     ylim = c(min(c(only_warm, only_evol))-0.0001, max(c(only_warm, only_evol))+0.0001))
abline(h = c(-0.0003, 0, 0.0003, 0.0006))
legend("topleft",
       legend = paste0("r = ", round(cor(dat$only_warm, dat$damage_present), 2)),
       bty = "n")

plot(only_evol ~ damage_present, dat,
     pch = 21, cex = 1, bg = "grey", lwd = 0.25,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     xlim = c(min(damage_present)-0.0001, max(damage_present)+0.0001),
     ylim = c(min(c(only_warm, only_evol))-0.0001, max(c(only_warm, only_evol))+0.0001))
abline(h = c(-0.0003, 0, 0.0003, 0.0006))
legend("topleft",
       legend = paste0("r = ", round(cor(dat$only_evol, dat$damage_present), 2)),
       bty = "n")
dev.off()
