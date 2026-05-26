#### Setup ----

source("packages_functions_data.R")

#### Priors and starting values ----

dev_rate_prior <- prior(normal(45, 4), nlpar = "Tmax") +
  prior(normal(10, 4), nlpar = "Tmin") +
  prior(normal(35, 4), nlpar = "Topt") +
  prior(normal(0.06, 0.06), nlpar = "Ropt", lb = 0)

init <- parse(text=paste0("list(b_Tmax = as.array(rnorm(3, 45, 0.1)),
  b_Tmin = as.array(rnorm(3, 10, 0.1)),
  b_Topt = as.array(rnorm(3, 35,0.1)),
  b_Ropt = as.array(rnorm(3, 0.06, 0.01)))"))
set.seed(6653);inits_list <- list(eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init))

#### Fit model ----

dev_rate_model_global <- brm(bf(dev.rate ~ log(LRF(temperature, Tmin, Tmax, Topt, Ropt * exp(randomRopt))),
                      Tmin + Tmax + Topt + Ropt ~ 0 + selection.regime,
                      randomRopt ~ 0 + (1 | replicate) + (1 | replicate:temperature.factor),
                         nl = TRUE),
                      data = dat,
                      prior = dev_rate_prior,
                      stanvars = stanvars,
                      family = lognormal(),
                      init = inits_list,
                      control = list(adapt_delta = 0.99, max_treedepth = 15),
                      warmup = 2000,
                      iter = 5000,
                      thin = 10,
                      seed = 1272,
                      cores = 8,
                      chains = 8,
                      file = "../models/dev_rate_model_global")

plot(dev_rate_model_global)
summary(dev_rate_model_global)

#### Posterior predictive check ----

pdf("../figures/raw/fig_S7_ppc_dev_rate_global.pdf", height = 10, width = 10, pointsize = 3)
set.seed(9536); print(pp_check(dev_rate_model_global, type = "ecdf_overlay_grouped", group = "replicate"))
dev.off()

#### Extract posteriors ----

dev_rate_model_global_posteriors <- as.data.frame(dev_rate_model_global)

#### Plot TPCs ----

pdf("../figures/raw/fig_5_dev_rate_TPC_global.pdf", height = 2, width = 2, pointsize = 3); set.seed(5285);  {

  plot(NA,
       ylab = "Development rate", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 0.065), xlim = c(10, 50),
       main = "Global", cex.main = 1.75)

  for(i in unique(dat$selection.regime)){

    tempdat <- subset(dat, selection.regime == i & both.sexes == 1)

    points(dev.rate ~ temperature,
           data = aggregate(dev.rate ~ temperature, tempdat, "geomean"),
           pch = 21, cex = 1.5,
           lwd = 0.5,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])

    TPC <- parse(text = paste0("curve(LRF(temp = x,
    Tmin = point_estimate(dev_rate_model_global_posteriors$`b_Tmin_selection.regime", i, "`, centrality ='MAP'),
    Tmax = point_estimate(dev_rate_model_global_posteriors$`b_Tmax_selection.regime", i, "`, centrality ='MAP'),
    Topt = point_estimate(dev_rate_model_global_posteriors$`b_Topt_selection.regime", i, "`, centrality ='MAP'),
      Ropt = point_estimate(dev_rate_model_global_posteriors$`b_Ropt_selection.regime", i, "`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = 1)"))

    eval(TPC)
  }

  ## y-axis
  ytick <- c(0, 0.025, 0.05)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)

}; dev.off()

#### Plot posterior estimates ----

pdf("../figures/raw/dev_rate_global_posteriors.pdf", height = 2.5, width = 1.5, pointsize = 3); for(parameter in c("Tmin", "Topt", "Tmax", "Ropt")) {
  iter <- 1

  ## Vectors for x- and y-lim values
  minima <- c()
  maxima <- c()
  for(i in c("anc", "cold", "hot")){
      eval(parse(text = paste0("minima <- c(minima, as.numeric(ci(dev_rate_model_global_posteriors$`b_", parameter, "_selection.regime", i, "`, type ='HDI', 0.95)[2]))
                                                                    maxima <- c(maxima, as.numeric(ci(dev_rate_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, type ='HDI', 0.95)[3]))")))
  }

  ypos <- iter
  mean_estimate <- c()

  plot(NA,
       ylab = "", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(43, 51), xlim = c(min(minima), max(maxima)),
       main = parameter, cex.main = 1.75)

  for(i in c("anc", "cold", "hot")){

        eval(parse(text = paste0("{
        dens <- as.data.frame(density(dev_rate_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, n = 10000))
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

        MAP <- eval(parse(text = paste0("point_estimate(dev_rate_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, centrality ='MAP')")))
        hdpi <- eval(parse(text = paste0("ci(dev_rate_model_global_posteriors$`b_", parameter,"_selection.regime", i, "`, type ='HDI', 0.85)")))

        # HPDI interval
        arrows(y0 = 50-ypos, y1 = 50-ypos,
               x0 = as.numeric(hdpi[2]), x1 = as.numeric(hdpi[3]),
               lwd = 2, code = 0,
               col = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), "black"))(6)[4])

        # MAP point estimate
        points(MAP, y = 50-ypos,
               pch = 21, cex = 1.5,
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

#### Extract parameter estimates ----

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
trait <- "development_rate"

## Main effects
for (parameter in parameters) {
    for (selection_regime in unique(dat$selection.regime)) {

      effect <- "main"
      param_name <- paste0("b_", parameter, "_selection.regime", selection_regime)

      if (param_name %in% names(dev_rate_model_global_posteriors)) {
        samples <- dev_rate_model_global_posteriors[[param_name]]
        MAP <- point_estimate(samples, centrality='MAP')
        HPDI <- ci(samples, method='HDI', ci=0.9)

        posterior_estimates <- rbind(posterior_estimates,
                                     data.frame(trait=trait,
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

## Random effects
for(random_params in c("sd_replicate__randomRopt_Intercept", "sd_replicate:temperature.factor__randomRopt_Intercept", "sigma")){
  samples <- dev_rate_model_global_posteriors[[random_params]]
  MAP <- point_estimate(samples, centrality='MAP')
  HPDI <- ci(samples, method='HDI', ci=0.9)

  parameter <- ifelse(random_params == "sd_replicate__randomRopt_Intercept", "replicate_sd",
                      ifelse(random_params == "sd_replicate:temperature.factor__randomRopt_Intercept", "replicateXtemp_sd", "resid_sd"))
  posterior_estimates <- rbind(posterior_estimates,
                               data.frame(trait=trait,
                                          parameter=parameter,
                                          effect="random",
                                          selection.regime=selection_regime,
                                          MAP=MAP,
                                          hpdi.lower=HPDI$CI_low,
                                          hpdi.upper=HPDI$CI_high,
                                          stringsAsFactors=FALSE))
}

## Save estimates
write.table(posterior_estimates, file = "../output/parameter_estimates_dev_rate_global.txt", sep = "\t", row.names = FALSE, quote = F)
