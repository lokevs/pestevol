#### Load packages, functions, and data ----
source("packages_functions_data.R")

#### Informative priors from development rate model ----

## Import Tmin/Tmax estimates from development rate
parameter_estimates_dev_rate <- subset(
  read.delim("../output/parameter_estimates_dev_rate_global.txt"),
  effect == "main" & parameter == "Tmax" | effect == "main" & parameter == "Tmin"
)

## Build parameter-specific priors from development rate posteriors
specific_priors_list <- list()
for (i in 1:nrow(parameter_estimates_dev_rate)) {
  nlpar_name <- parameter_estimates_dev_rate$parameter[i]
  selection_val <- parameter_estimates_dev_rate$selection.regime[i]
  map_estimate <- parameter_estimates_dev_rate$MAP[i]
  prior_sd <- 1
  coef_name <- paste0("b_", nlpar_name, "_selection.regime", selection_val)

  current_prior <- eval(parse(text = paste0("prior(normal(", map_estimate, ", ", prior_sd, "),
                         nlpar = ", nlpar_name, ",
                         coef = ", gsub("b_.*?_", "", coef_name), ")")))

  specific_priors_list[[length(specific_priors_list) + 1]] <- current_prior
}

specific_priors_list <- Reduce(`+`, specific_priors_list)

## Combine with remaining priors
growth_rate_prior <- prior(normal(35, 3), nlpar = "Topt") +
  prior(normal(0.00008, 0.00001), nlpar = "Ropt", lb = 0) +
  specific_priors_list

#### Fit growth rate model ----

## Starting values (randomized)
init <- parse(text=paste0("list(b_Tmax = as.array(rnorm(3, 45, 0.1)),
  b_Tmin = as.array(rnorm(3, 10, 0.1)),
  b_Topt = as.array(rnorm(3, 35,0.1)),
  b_Tmax = as.array(rnorm(3, 45,0.1)),
  b_Ropt = as.array(rnorm(3, 0.06, 0.01)))"))
set.seed(6863);inits_list <- list(eval(init), eval(init), eval(init), eval(init))

## Model
growth_rate_model_global <- brm(bf(growth.rate ~ log(LRF(temperature, Tmin, Tmax, Topt, Ropt * exp(randomRopt))),
                            Tmin  + Topt  + Tmax  + Ropt ~ 0 + selection.regime,
                            randomRopt ~ 0 + (1 | replicate) + (1 | replicate:temperature.factor),
                            nl = TRUE),
                         data = dat,
                         prior = growth_rate_prior,
                         stanvars = stanvars,
                         family = lognormal(),
                         init = inits_list,
                         control = list(adapt_delta = 0.9, max_treedepth = 15),
                         warmup = 2000,
                         iter = 5000,
                         thin = 10,
                         seed = 3567,
                         cores = 4,
                         chains = 4,
                         file = "../models/growth_rate_model_global")

plot(growth_rate_model_global)
summary(growth_rate_model_global)

#### Posterior predictive check ----
pdf("../figures/raw/fig_S9_ppc_growth_rate_global.pdf", height = 10, width = 10, pointsize = 3); set.seed(2491); pp_check(growth_rate_model_global, type = "ecdf_overlay_grouped", group = "replicate"); dev.off()

#### Plot growth rate TPCs ----

growth_rate_model_global_posteriors <- as.data.frame(growth_rate_model_global)

pdf("../figures/raw/fig_5_growth_rate_TPCs_global.pdf", height = 2, width = 2, pointsize = 3); set.seed(5285);{

plot(NA,
     ylab = "Growth rate", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
     ylim = c(0, 0.00015), xlim = c(10, 50),
     main = i, cex.main = 1.75)

for(i in c("anc", "cold", "hot")) {

    tempdat <- subset(dat, selection.regime == i & both.sexes == 1)

    points(growth.rate ~ temperature,
           data = aggregate(growth.rate ~ temperature, tempdat, "geomean"),
           pch = 21, cex = 1.5,
           lwd = 0.5,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])

    TPC <- parse(text = paste0("curve(LRF(temp = x,
    Tmin = point_estimate(growth_rate_model_global_posteriors$`b_Tmin_selection.regime", i, "`, centrality ='MAP'),
    Tmax = point_estimate(growth_rate_model_global_posteriors$`b_Tmax_selection.regime", i, "`, centrality ='MAP'),
    Topt = point_estimate(growth_rate_model_global_posteriors$`b_Topt_selection.regime", i, "`, centrality ='MAP'),
      Ropt = point_estimate(growth_rate_model_global_posteriors$`b_Ropt_selection.regime", i, "`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = 1)"))

    eval(TPC)

  ytick <- c(0, 0.00005, 0.0001, 0.00015)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick * 1000, pos = 2, xpd = TRUE, cex=1.5)

}}; dev.off()

#### Plot posterior distributions ----
pdf("../figures/raw/growth_rate_posteriors_global.pdf", height = 2.5, width = 1.5, pointsize = 3); for(parameter in c("Tmin", "Topt", "Tmax", "Ropt")) {
  iter <- 1

  minima <- c()
  maxima <- c()
  for(i in c("anc", "cold", "hot")){
      eval(parse(text = paste0("minima <- c(minima, as.numeric(ci(growth_rate_model_global_posteriors$`b_", parameter, "_selection.regime", i, "`, type ='HDI', 0.95)[2]))
                                                                    maxima <- c(maxima, as.numeric(ci(growth_rate_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, type ='HDI', 0.95)[3]))")))
    }

  ypos <- iter
  mean_estimate <- c()

  plot(NA,
       ylab = "", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(43, 51), xlim = c(min(minima), max(maxima)),
       main = parameter, cex.main = 1.75)

  for(i in c("anc", "cold", "hot")){

      eval(parse(text = paste0("{
        dens <- as.data.frame(density(growth_rate_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, n = 10000))
        sorted_density <- sort(dens$y, decreasing = TRUE)
        cumulative_density <- cumsum(sorted_density) / sum(dens$y)
        cutoff <- sorted_density[which(cumulative_density >= 0.90)[1]]
        dens_trimmed <- dens
        dens_trimmed$y <- ifelse(dens$y >= cutoff, dens$y, NA)
        dens_trimmed <- dens_trimmed[!is.na(dens_trimmed$y),]

        polygon(x = c(min(dens_trimmed$x), dens_trimmed$x, max(dens_trimmed$x)),
     y = (c(0, dens_trimmed$y, 0)) / max(dens_trimmed$y) + 50 - ypos,
                lwd = 0.25,
                border =  colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), 'black'))(6)[4],
                col = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), 'white'))(6)[4])
                                 }")))

      MAP <- eval(parse(text = paste0("point_estimate(growth_rate_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, centrality ='MAP')")))
      hdpi <- eval(parse(text = paste0("ci(growth_rate_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, type ='HDI', 0.85)")))

      arrows(y0 = 50-ypos, y1 = 50-ypos,
             x0 = as.numeric(hdpi[2]), x1 = as.numeric(hdpi[3]),
             lwd = 2, code = 0,
             col = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), "black"))(6)[4])

      points(MAP, y = 50-ypos,
             pch =   21,
             cex = 1.5,
             lwd = 0.5,
             bg = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), "white"))(6)[1],
             col = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), 'black'))(6)[4])

      ypos = ypos + 1
      mean_estimate <- c(mean_estimate, as.numeric(MAP))

    ypos = ypos + 1.5
    iter = iter + 1

  }

  abline(v = mean(mean_estimate), lty = 3)

}; dev.off()

#### Extract and export posterior estimates ----

posterior_estimates <- data.frame(trait=character(),
                                  parameter=character(),
                                  effect=character(),
                                  origin=character(),
                                  selection.regime=character(),
                                  replicate=character(),
                                  MAP=numeric(),
                                  hpdi.lower=numeric(),
                                  hpdi.upper=numeric(),
                                  stringsAsFactors=FALSE)

parameters <- c("Tmin", "Topt", "Tmax", "Ropt")

## Main effects
for (parameter in parameters) {
  for (selection_regime in unique(dat$selection.regime)) {

    effect <- "main"
    param_name <- paste0("b_", parameter, "_selection.regime", selection_regime)

    if (param_name %in% names(growth_rate_model_global_posteriors)) {
      samples <- growth_rate_model_global_posteriors[[param_name]]
      MAP <- point_estimate(samples, centrality='MAP')
      HPDI <- ci(samples, method='HDI', ci=0.9)

      posterior_estimates <- rbind(posterior_estimates,
                                   data.frame(trait="growth_rate",
                                              parameter=parameter,
                                              effect=effect,
                                              selection.regime=selection_regime,
                                              MAP=MAP,
                                              hpdi.lower=HPDI$CI_low,
                                              hpdi.upper=HPDI$CI_high,
                                              stringsAsFactors=FALSE))
    } else {
      message(paste("Parameter", param_name, "not found in posterior samples."))
    }
  }
}

## Random effects and residual SD
for(random_params in c("sd_replicate__randomRopt_Intercept", "sd_replicate:temperature.factor__randomRopt_Intercept", "sigma")){
  samples <- growth_rate_model_global_posteriors[[random_params]]
  MAP <- point_estimate(samples, centrality='MAP')
  HPDI <- ci(samples, method='HDI', ci=0.9)

  parameter <- ifelse(random_params == "sd_replicate__randomRopt_Intercept", "replicate_sd",
                      ifelse(random_params == "sd_replicate:temperature.factor__randomRopt_Intercept", "replicateXtemp_sd", "resid_sd"))
  posterior_estimates <- rbind(posterior_estimates,
                               data.frame(trait="growth_rate",
                                          parameter=parameter,
                                          effect="random",
                                          selection.regime=selection_regime,
                                          MAP=MAP,
                                          hpdi.lower=HPDI$CI_low,
                                          hpdi.upper=HPDI$CI_high,
                                          stringsAsFactors=FALSE))
}

write.table(posterior_estimates, file = "../output/parameter_estimates_growth_rate_global.txt", sep = "\t", row.names = FALSE, quote = F)
