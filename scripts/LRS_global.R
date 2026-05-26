#### Load requirements ----
source("packages_functions_data.R")

#### Specify priors and starting values ----
fecundity_model_prior <- prior(normal(29, 5), nlpar="Topt", lb = 0) +
  prior(normal(70,25), nlpar="fecmax", lb = 0) +
  prior(lognormal(log(22), log(1.1)), nlpar="breadth", lb = 0) +
  prior(lognormal((2.5),(2.5)), nlpar="fatness", lb = 0) +
  prior(lognormal(log(0.03), log(1.2)), nlpar="a", lb = 0) +
  prior(normal(28, 2), nlpar="h") +
  prior(normal(-4.5, 2), nlpar="k")+
  prior(normal(0, 2.5), dpar="shape")

## Starting values (randomized)
init <- parse(text=paste0("list(b_Topt = as.array(rnorm(3, 29, 1.5)),
  b_fecmax = as.array(rnorm(3, 70, 1.5)),
  b_breadth = as.array(rlnorm(3, log(25), log(1.1))),
  b_fatness = as.array(rlnorm(3, log(4), log(1.5))))"))
set.seed(1646);inits_list <- list(eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init))

#### Fit model ----
offspring_model_global <- brm(bf(adult.offspring ~ log(fec_function(temperature, Topt, fecmax, breadth, fatness)) + ranefffec,
                          nlf(zi ~ a * (temperature - h)^2 + k  + ranefffert),
                          shape ~ 0 + temperature.factor:selection.regime,
                          Topt + fecmax + breadth + a + h + k + fatness ~ 0 + selection.regime,
                          ranefffert + ranefffec ~ 0 + (1 | replicate) +  (1 | replicate:temperature.factor),
                          nl = T),
                       data = dat,
                       family = zero_inflated_negbinomial,
                       prior = fecundity_model_prior,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       warmup = 2000,
                       init = inits_list,
                       stanvars = stanvars,
                       iter = 5000,
                       thin = 1,
                       seed = 1688,
                       cores = 6,
                       chains = 6,
                       file = "../models/offspring_model_global")

plot(offspring_model_global)
summary(offspring_model_global)

#### Posterior predictive check ----
pdf("../figures/raw/fig_S8_ppc_LRS_global.pdf", height = 10, width = 10, pointsize = 3); set.seed(7434); print(pp_check(offspring_model_global, type = "ecdf_overlay_grouped", group = "replicate")); dev.off()

#### Extract posteriors ----
offspring_model_global_posteriors <- as.data.frame(offspring_model_global)

#### Plot fecundity and non-viability curves ----
pdf("../figures/raw/fig_5_fecundity_curves_global.pdf", height = 2, width = 2, pointsize = 3); set.seed(1589); {

plot(NA,
     ylab = "Fecundity", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
     ylim = c(0, 120), xlim = c(10, 50),
     main = i, cex.main = 1.75)

for(i in c("anc", "cold", "hot")) {

    tempdat <- subset(dat, selection.regime == i & both.sexes == 1)

    # Raw data
    points(adult.offspring ~ temperature,
           data = aggregate(adult.offspring ~ temperature, tempdat[tempdat$adult.offspring != 0,], "geomean"),
           pch = 21, cex = 1.5,
           lwd = 0.5,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])

    fec <- parse(text = paste0("curve(fec_function(temp = x,
    Topt = point_estimate(offspring_model_global_posteriors$`b_Topt_selection.regime", i, "`, centrality ='MAP'),
    fecmax = point_estimate(offspring_model_global_posteriors$`b_fecmax_selection.regime", i, "`, centrality ='MAP'),
    breadth = point_estimate(offspring_model_global_posteriors$`b_breadth_selection.regime", i, "`, centrality ='MAP'),
      fatness = point_estimate(offspring_model_global_posteriors$`b_fatness_selection.regime", i, "`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = 1)"))

    eval(fec)

    viability <- parse(text = paste0("curve(invlogit(viab_function(temp = x,
    a = point_estimate(offspring_model_global_posteriors$`b_a_selection.regime", i, "`, centrality ='MAP'),
    h = point_estimate(offspring_model_global_posteriors$`b_h_selection.regime", i, "`, centrality ='MAP'),
    k = point_estimate(offspring_model_global_posteriors$`b_k_selection.regime", i, "`, centrality ='MAP'))) * 120,
      10, 50, add = T, col = colorRampPalette(c(unique(tempdat$colour), 'black'))(6)[3], lty = 1, lwd = 1)"))

    eval(viability)

  ytick <- c(0, 60, 120)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)

  }
  }; dev.off()

#### Plot posterior estimates ----
pdf("../figures/raw/fecundity_posteriors_global.pdf", height = 2.5, width = 1.5, pointsize = 3); for(parameter in c("Topt", "fecmax", "breadth", "fatness", "h", "k", "a")){

      minima <- c()
      maxima <- c()
      for(i in c("anc", "cold", "hot")){
          eval(parse(text = paste0("minima <- c(minima, as.numeric(ci(offspring_model_global_posteriors$`b_", parameter, "_selection.regime", i, "`, type ='HDI', 0.95)[2]))
                                                                    maxima <- c(maxima, as.numeric(ci(offspring_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, type ='HDI', 0.95)[3]))")))
        }

  iter <- 1
  ypos <- iter
  mean_estimate <- c()

  plot(NA,
       ylab = "", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(43, 51), xlim = c(min(minima), max(maxima)),
       main = parameter, cex.main = 1.75)

  for(i in c("anc", "cold", "hot")){

      eval(parse(text = paste0("{
        dens <- as.data.frame(density(offspring_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, n = 10000))
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

      MAP <- eval(parse(text = paste0("point_estimate(na.omit(offspring_model_global_posteriors$`b_", parameter, "_selection.regime", i, "`), centrality ='MAP')")))
      hdpi <- eval(parse(text = paste0("ci(na.omit(offspring_model_global_posteriors$`b_", parameter, "_selection.regime", i, "`), type ='HDI', 0.85)")))

      arrows(y0 = 50-ypos, y1 = 50-ypos,
             x0 = as.numeric(hdpi[2]), x1 = as.numeric(hdpi[3]),
             lwd = 2, code = 0,
             col = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), "black"))(6)[4])

      points(MAP, y = 50-ypos,
             pch = 21,
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

#### Extract and save parameter estimates ----
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

## Fecundity parameters
fecundity_parameters <- c("Topt", "fecmax", "breadth", "fatness")

for (parameter in fecundity_parameters) {
    for (selection_regime in unique(dat$selection.regime)) {

      effect <- "main"
      param_name <- paste0("b_", parameter, "_selection.regime", selection_regime)

      if (param_name %in% names(offspring_model_global_posteriors)) {
        samples <- offspring_model_global_posteriors[[param_name]]
        MAP <- point_estimate(samples, centrality='MAP')
        HPDI <- ci(samples, method='HDI', ci=0.9)

        posterior_estimates <- rbind(posterior_estimates,
                                     data.frame(trait="fecundity",
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

## Fecundity random effects
for(random_params in c("sd_replicate__ranefffec_Intercept", "sd_replicate:temperature.factor__ranefffec_Intercept")){
  samples <- offspring_model_global_posteriors[[random_params]]
  MAP <- point_estimate(samples, centrality='MAP')
  HPDI <- ci(samples, method='HDI', ci=0.9)

  parameter <- ifelse(random_params == "sd_replicate__ranefffec_Intercept", "fec_replicate_sd", "fec_replicateXtempresid_sd")
  posterior_estimates <- rbind(posterior_estimates,
                               data.frame(trait="fecundity",
                                          parameter=parameter,
                                          effect="random",
                                          selection.regime=selection_regime,
                                          MAP=MAP,
                                          hpdi.lower=HPDI$CI_low,
                                          hpdi.upper=HPDI$CI_high,
                                          stringsAsFactors=FALSE))
}

## Viability parameters
viability_parameters <- c("h", "k", "a")

for (parameter in viability_parameters) {
  for (selection_regime in unique(dat$selection.regime)) {
    effect <- "main"
    param_name <- paste0("b_", parameter, "_selection.regime", selection_regime)

    if (param_name %in% names(offspring_model_global_posteriors)) {
      samples <- offspring_model_global_posteriors[[param_name]]
      MAP <- point_estimate(samples, centrality='MAP')
      HPDI <- ci(samples, method='HDI', ci=0.9)

      posterior_estimates <- rbind(posterior_estimates,
                                   data.frame(trait="viability",
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

## Viability random effects
for(random_params in c("sd_replicate__ranefffert_Intercept", "sd_replicate:temperature.factor__ranefffert_Intercept")){
  samples <- offspring_model_global_posteriors[[random_params]]
  MAP <- point_estimate(samples, centrality='MAP')
  HPDI <- ci(samples, method='HDI', ci=0.9)

  parameter <- ifelse(random_params == "sd_replicate__ranefffert_Intercept", "fert_replicate_sd", "fert_replicateXtemp_sd")

  posterior_estimates <- rbind(posterior_estimates,
                               data.frame(trait="fecundity",
                                          parameter=parameter,
                                          effect="random",
                                          selection.regime=selection_regime,
                                          MAP=MAP,
                                          hpdi.lower=HPDI$CI_low,
                                          hpdi.upper=HPDI$CI_high,
                                          stringsAsFactors=FALSE))
}

write.table(posterior_estimates, file = "../output/parameter_estimates_viab_fec_global.txt", sep = "\t", row.names = FALSE, quote = F)
