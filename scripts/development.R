#### Load requirements ----
source("packages_functions_data.R")

#### Priors and starting values ----

dev_rate_prior <- prior(normal(45, 4), nlpar = "Tmax") +
  prior(normal(10, 4), nlpar = "Tmin") +
  prior(normal(35, 4), nlpar = "Topt") +
  prior(normal(0.06, 0.06), nlpar = "Ropt", lb = 0)

## Randomized starting values
init <- parse(text=paste0("list(b_Tmax = as.array(rnorm(9, 45, 0.1)),
  b_Tmin = as.array(rnorm(9, 10, 0.1)),
  b_Topt = as.array(rnorm(9, 35,0.1)),
  b_Ropt = as.array(rnorm(9, 0.06, 0.01)))"))
set.seed(1705);inits_list <- list(eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init))

#### Fit model ----

dev_rate_model <- brm(bf(dev.rate ~ log(LRF(temperature, Tmin, Tmax, Topt, Ropt * exp(randomRopt))),
                      Tmin + Tmax + Topt ~ 0 + origin:selection.regime + random.dummy:(1 | replicate),
                      Ropt ~ 0 + origin:selection.regime,
                      randomRopt ~ 0 + random.dummy:(1 | replicate),
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
                      seed = 5671,
                      cores = 8,
                      chains = 8,
                      file = "../models/dev_rate_model")

plot(dev_rate_model)
summary(dev_rate_model)

#### Posterior predictive check ----

pdf("../figures/raw/fig_S4_ppc_dev_rate.pdf", height = 10, width = 10, pointsize = 3)
set.seed(9536); print(pp_check(dev_rate_model, type = "ecdf_overlay_grouped", group = "replicate"))
dev.off()

#### Extract posteriors ----

dev_rate_model_posteriors <- as.data.frame(dev_rate_model)

#### Plot TPCs with all data ----

pdf("../figures/raw/fig_S1_dev_rate_TPCs_all_data.pdf", height = 2, width = 2, pointsize = 3); set.seed(5285); for(i in c("anc", "cold", "hot")) {

  plot(NA,
       ylab = "Development rate", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 0.065), xlim = c(10, 50),
       main = i, cex.main = 1.75)

  for(g in unique(subset(dat, selection.regime == i)$origin)){

    tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

    points(jitter(dev.rate) ~ jitter(temperature),
           data = tempdat,
           pch = unique(tempdat$point), cex = 1,
           lwd = 0.5,
           bg = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[4],
           col = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[2])

    if(unique(tempdat$origin) == "bra"){linetype <- 1}
    if(unique(tempdat$origin) == "yem"){linetype <- 2}
    if(unique(tempdat$origin) == "ca"){linetype <- 3}

    TPC <- parse(text = paste0("curve(LRF(temp = x,
    Tmin = point_estimate(dev_rate_model_posteriors$`b_Tmin_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    Tmax = point_estimate(dev_rate_model_posteriors$`b_Tmax_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    Topt = point_estimate(dev_rate_model_posteriors$`b_Topt_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
      Ropt = point_estimate(dev_rate_model_posteriors$`b_Ropt_origin", g, ":selection.regime", i, "`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = linetype)"))

    eval(TPC)
  }

  ## Group means
  groupmeans_regime <- aggregate(dev.rate ~ replicate + temperature + origin + colour + point, data = subset(dat, selection.regime == i & both.sexes == 1 & adult.offspring != 0), "geomean")
  groupmeans_regime <- aggregate(dev.rate ~ temperature + origin + colour + point, data = groupmeans_regime, "geomean")

  for(g in unique(groupmeans_regime$origin)){

    tempdat <- subset(groupmeans_regime, origin == g)

    points(dev.rate ~ temperature,
           data = tempdat,
           pch = unique(tempdat$point), cex = 2,
           lwd = 0.5,
           lty = 3,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  ## Y-axis
  ytick <- c(0, 0.025, 0.05)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)

}; legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n"); dev.off()

#### Plot TPCs clean ----

pdf("../figures/raw/fig_S3_dev_rate_TPCs.pdf", height = 2, width = 2, pointsize = 3); set.seed(5285); for(i in c("anc", "cold", "hot")) {

  plot(NA,
       ylab = "Development rate", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 0.065), xlim = c(10, 50),
       main = i, cex.main = 1.75)

  abline(h = c(0, 0.025, 0.05), col = "grey40", lwd = 0.5, lty = 3)

  for(g in unique(subset(dat, selection.regime == i)$origin)){

    tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

    if(unique(tempdat$origin) == "bra"){linetype <- 1}
    if(unique(tempdat$origin) == "yem"){linetype <- 2}
    if(unique(tempdat$origin) == "ca"){linetype <- 3}

    TPC <- parse(text = paste0("curve(LRF(temp = x,
    Tmin = point_estimate(dev_rate_model_posteriors$`b_Tmin_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    Tmax = point_estimate(dev_rate_model_posteriors$`b_Tmax_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    Topt = point_estimate(dev_rate_model_posteriors$`b_Topt_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
      Ropt = point_estimate(dev_rate_model_posteriors$`b_Ropt_origin", g, ":selection.regime", i, "`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = linetype)"))

    eval(TPC)
  }

  ## Group means
  groupmeans_regime <- aggregate(dev.rate ~ replicate + temperature + origin + colour + point, data = subset(dat, selection.regime == i & both.sexes == 1 & adult.offspring != 0), "geomean")
  groupmeans_regime <- aggregate(dev.rate ~ temperature + origin + colour + point, data = groupmeans_regime, "geomean")

  for(g in unique(groupmeans_regime$origin)){

    tempdat <- subset(groupmeans_regime, origin == g)

    points(dev.rate ~ temperature,
           data = tempdat,
           pch = unique(tempdat$point), cex = 2,
           lwd = 0.5,
           lty = 3,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  ## Y-axis
  ytick <- c(0, 0.025, 0.05)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)

}; legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n"); dev.off()

#### Plot posterior estimates ----

pdf("../figures/raw/fig_S2_dev_rate_posteriors.pdf", height = 2.5, width = 1.5, pointsize = 3); for(parameter in c("Tmin", "Topt", "Tmax", "Ropt")) {
  iter <- 1

  ## Compute axis limits from 95% HDIs
  minima <- c()
  maxima <- c()
  for(i in c("anc", "cold", "hot")){
    for(g in unique(subset(dat, selection.regime == i)$origin)){
      eval(parse(text = paste0("minima <- c(minima, as.numeric(ci(dev_rate_model_posteriors$`b_", parameter, "_origin", g, ":selection.regime", i, "`, type ='HDI', 0.95)[2]))
                                                                    maxima <- c(maxima, as.numeric(ci(dev_rate_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, type ='HDI', 0.95)[3]))")))
    }
  }

  ypos <- iter
  mean_estimate <- c()

  plot(NA,
       ylab = "", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(37, 51), xlim = c(min(minima), max(maxima)),
       main = parameter, cex.main = 1.75)

  for(i in c("anc", "cold", "hot")){
    for(g in unique(subset(dat, selection.regime == i)$origin)){

        eval(parse(text = paste0("{
        dens <- as.data.frame(density(dev_rate_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, n = 10000))
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

        MAP <- eval(parse(text = paste0("point_estimate(dev_rate_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, centrality ='MAP')")))
        hdpi <- eval(parse(text = paste0("ci(dev_rate_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, type ='HDI', 0.85)")))

        # 85% HPDI
        arrows(y0 = 50-ypos, y1 = 50-ypos,
               x0 = as.numeric(hdpi[2]), x1 = as.numeric(hdpi[3]),
               lwd = 2, code = 0,
               col = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), "black"))(6)[4])

        # MAP point estimate
        points(MAP, y = 50-ypos,
               pch = unique(subset(dat, origin == g)$point), cex = 1.5,
               lwd = 0.5,
               bg = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), "white"))(6)[1],
               col = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), 'black'))(6)[4])

        ypos = ypos + 1
        mean_estimate <- c(mean_estimate, as.numeric(MAP))

    }

    ypos = ypos + 1.5
    iter = iter + 1

  }

  abline(v = mean(mean_estimate), lty = 3)

  }; dev.off()

#### Extract posterior estimates and uncertainties ----

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
  for (origin in unique(dat$origin)) {
    for (selection_regime in unique(dat$selection.regime)) {

      effect <- "main"
      param_name <- paste0("b_", parameter, "_origin", origin, ":selection.regime", selection_regime)

      if (param_name %in% names(dev_rate_model_posteriors)) {
        samples <- dev_rate_model_posteriors[[param_name]]
        MAP <- point_estimate(samples, centrality='MAP')
        HPDI <- ci(samples, method='HDI', ci=0.9)

        posterior_estimates <- rbind(posterior_estimates,
                                     data.frame(trait=trait,
                                                parameter=parameter,
                                                effect=effect,
                                                origin=origin,
                                                selection.regime=selection_regime,
                                                replicate=NA,
                                                MAP=MAP,
                                                hpdi.lower=HPDI$CI_low,
                                                hpdi.upper=HPDI$CI_high,
                                                stringsAsFactors=FALSE))
      } else {
        message(paste("Parameter", param_name, "not found in posterior samples."))
      }
    }
  }
}

## Random effects (Tmin, Topt, Tmax)
for (parameter in c("Tmin", "Topt", "Tmax")) {
  for (replicate in unique(dat$replicate)) {

    random_dummy <- unique(dat$random.dummy[dat$replicate == replicate])

    if (random_dummy == 1) {

      effect <- "random"
      origin <- unique(dat$origin[dat$replicate == replicate])
      selection_regime <- unique(dat$selection.regime[dat$replicate == replicate])

      fixed_param_name <- paste0("b_", parameter, "_origin", origin, ":selection.regime", selection_regime)
      random_param_name <- paste0("r_replicate__", parameter, "[", replicate, ",Intercept]")

      if (fixed_param_name %in% names(dev_rate_model_posteriors) && random_param_name %in% names(dev_rate_model_posteriors)) {

        total_samples <- dev_rate_model_posteriors[[fixed_param_name]] + dev_rate_model_posteriors[[random_param_name]]
        MAP <- point_estimate(total_samples, centrality='MAP')
        HPDI <- ci(total_samples, method='HDI', ci=0.9)

        posterior_estimates <- rbind(posterior_estimates,
                                     data.frame(trait=trait,
                                                parameter=parameter,
                                                effect=effect,
                                                origin=origin,
                                                selection.regime=selection_regime,
                                                replicate=replicate,
                                                MAP=MAP,
                                                hpdi.lower=HPDI$CI_low,
                                                hpdi.upper=HPDI$CI_high,
                                                stringsAsFactors=FALSE))
      } else {
        message(paste("Parameters", fixed_param_name, "or", random_param_name, "not found in posterior samples."))
      }
    }
  }
}

## Random effects (Ropt via exp(randomRopt))
parameter <- "randomRopt"
for (replicate in unique(dat$replicate)) {

  random_dummy <- unique(dat$random.dummy[dat$replicate == replicate])

  if (random_dummy == 1) {

    effect <- "random"
    origin <- unique(dat$origin[dat$replicate == replicate])
    selection_regime <- unique(dat$selection.regime[dat$replicate == replicate])

    fixed_param_name <- paste0("b_Ropt_origin", origin, ":selection.regime", selection_regime)
    random_param_name <- paste0("r_replicate__", parameter, "[", replicate, ",Intercept]")

    if (fixed_param_name %in% names(dev_rate_model_posteriors) && random_param_name %in% names(dev_rate_model_posteriors)) {

      # Ropt * exp(randomRopt)
      total_samples <- dev_rate_model_posteriors[[fixed_param_name]] * exp(dev_rate_model_posteriors[[random_param_name]])
      MAP <- point_estimate(total_samples, centrality='MAP')
      HPDI <- ci(total_samples, method='HDI', ci=0.9)

      posterior_estimates <- rbind(posterior_estimates,
                                   data.frame(trait=trait,
                                              parameter="Ropt",
                                              effect=effect,
                                              origin=origin,
                                              selection.regime=selection_regime,
                                              replicate=replicate,
                                              MAP=MAP,
                                              hpdi.lower=HPDI$CI_low,
                                              hpdi.upper=HPDI$CI_high,
                                              stringsAsFactors=FALSE))
    } else {
      message(paste("Parameters", fixed_param_name, "or", random_param_name, "not found in posterior samples."))
    }
  }
}

#### Save parameter estimates ----

posterior_estimates_dev_rate <- posterior_estimates
write.table(posterior_estimates_dev_rate, file = "../output/parameter_estimates_dev_rate.txt", sep = "\t", row.names = FALSE, quote = F)
