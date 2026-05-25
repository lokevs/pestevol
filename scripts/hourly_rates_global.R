#### Setup ----
source("packages_functions_data.R")

#### Load global (regime-only) parameter estimates ----
TPC_parameters_global <- rbind(read.table("../output/parameter_estimates_viab_fec_global.txt", header = TRUE),
                               read.table("../output/parameter_estimates_dev_rate_global.txt", header = TRUE),
                               read.table("../output/parameter_estimates_growth_rate_global.txt", header = TRUE))
TPC_parameters_global <- subset(TPC_parameters_global, effect == "main")

get_param <- function(tr, par, reg){
  TPC_parameters_global$MAP[TPC_parameters_global$trait == tr &
                              TPC_parameters_global$parameter == par &
                              TPC_parameters_global$selection.regime == reg]
}

#### Plot settings ----
regimes <- c("anc", "cold", "hot")
regime_colors <- c(anc = "darkgrey", cold = "dodgerblue", hot = "orangered")
temp_seq <- seq(5, 45, 0.1)

## Pre-compute curves for y-limits
dev_curves <- sapply(regimes, function(reg)
  LRF_hourly(temp_seq,
             Tmin = get_param("development_rate", "Tmin", reg),
             Tmax = get_param("development_rate", "Tmax", reg),
             Topt = get_param("development_rate", "Topt", reg),
             Ropt = get_param("development_rate", "Ropt", reg)))

growth_curves <- sapply(regimes, function(reg)
  LRF_hourly(temp_seq,
             Tmin = get_param("growth_rate", "Tmin", reg),
             Tmax = get_param("growth_rate", "Tmax", reg),
             Topt = get_param("growth_rate", "Topt", reg),
             Ropt = get_param("growth_rate", "Ropt", reg)))

viab_curves <- sapply(regimes, function(reg)
  viability_rate_hourly(temp_seq,
                        Tmin = get_param("development_rate", "Tmin", reg),
                        Tmax = get_param("development_rate", "Tmax", reg),
                        Topt = get_param("development_rate", "Topt", reg),
                        Ropt = get_param("development_rate", "Ropt", reg),
                        a = get_param("viability", "a", reg),
                        h = get_param("viability", "h", reg),
                        k = get_param("viability", "k", reg)))

ovi_curves <- sapply(regimes, function(reg){
  LRS_vals <- fec_function(temp_seq,
                           Topt = get_param("fecundity", "Topt", reg),
                           fecmax = get_param("fecundity", "fecmax", reg),
                           breadth = get_param("fecundity", "breadth", reg),
                           fatness = get_param("fecundity", "fatness", reg))
  if(reg == "anc") oviposition_rate_function_anc_hourly(temp_seq, LRS_vals)
  else             oviposition_rate_function_selection_hourly(temp_seq, LRS_vals)
})

#### 4-panel figure ----
pdf("../figures/raw/fig_5_hourly_rates_global.pdf", height = 1.5, width = 1.3, pointsize = 3)

## (1) Development rate
plot(NA, xlim = c(5, 45), ylim = c(0, max(dev_curves) * 1.05),
     xlab = "Temperature (°C)", ylab = "Hourly development rate",
     main = "Development", cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)
for(reg in regimes) lines(temp_seq, dev_curves[, reg], col = regime_colors[reg], lwd = 1)

## (2) Growth rate
plot(NA, xlim = c(5, 45), ylim = c(0, max(growth_curves) * 1.05),
     xlab = "Temperature (°C)", ylab = "Hourly growth rate",
     main = "Growth", cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)
for(reg in regimes) lines(temp_seq, growth_curves[, reg], col = regime_colors[reg], lwd = 1)

## (3) Viability rate
plot(NA, xlim = c(5, 45), ylim = c(0.99, 1),
     xlab = "Temperature (°C)", ylab = "Hourly viability probability",
     main = "Viability", cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)
for(reg in regimes) lines(temp_seq, viab_curves[, reg], col = regime_colors[reg], lwd = 1)

## (4) Oviposition rate
plot(NA, xlim = c(5, 45), ylim = c(0, max(ovi_curves) * 1.05),
     xlab = "Temperature (°C)", ylab = "Hourly oviposition rate",
     main = "Oviposition", cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)
for(reg in regimes) lines(temp_seq, ovi_curves[, reg], col = regime_colors[reg], lwd = 1)

legend("topleft", legend = c("ancestral", "cold", "hot"),
       col = regime_colors[regimes], lwd = 2, bty = "n", cex = 1.2)

dev.off()
