#### Load packages, functions, and data ----
source("packages_functions_data.R")

#### Informative priors from development rate model ----

parameter_estimates_dev_rate <- subset(read.delim("../output/parameter_estimates_dev_rate.txt"), effect == "main" & parameter == "Tmax" | effect == "main" & parameter == "Tmin")
specific_priors_list <- list()
for (i in 1:nrow(parameter_estimates_dev_rate)) {

  nlpar_name <- parameter_estimates_dev_rate$parameter[i]
  origin_val <- parameter_estimates_dev_rate$origin[i]
  selection_val <- parameter_estimates_dev_rate$selection.regime[i]
  map_estimate <- parameter_estimates_dev_rate$MAP[i]
  prior_sd <- 1
  coef_name <- paste0("b_", nlpar_name, "_origin", origin_val, ":selection.regime", selection_val)

  current_prior <- eval(parse(text = paste0("prior(normal(", map_estimate, ", ", prior_sd, "),
                         nlpar = ", nlpar_name, ",
                         coef = ", gsub("b_.*?_", "", coef_name), ")")))

  specific_priors_list[[length(specific_priors_list) + 1]] <- current_prior
}

specific_priors_list <- Reduce(`+`, specific_priors_list)

growth_rate_prior <- prior(normal(35, 3), nlpar = "Topt") +
  prior(normal(0.00008, 0.00001), nlpar = "Ropt", lb = 0) +
  prior(normal(0, 0.25), class="sd", group="replicate", nlpar="Tmax", lb = 0) +
  specific_priors_list

#### Fit growth rate model ----

## Starting values
init <- parse(text=paste0("list(b_Tmax = as.array(rnorm(9, 45, 0.1)),
  b_Tmin = as.array(rnorm(9, 10, 0.1)),
  b_Topt = as.array(rnorm(9, 35,0.1)),
  b_Ropt = as.array(rnorm(9, 0.06, 0.01)))"))
set.seed(6863);inits_list <- list(eval(init), eval(init), eval(init), eval(init))

## Model
growth_rate_model <- brm(bf(growth.rate ~ log(LRF(temperature, Tmin, Tmax, Topt, Ropt * exp(randomRopt))),
                            Tmin  + Topt  + Tmax ~ 0 + origin:selection.regime + random.dummy:(1 | replicate),
                            Ropt ~ 0 + origin:selection.regime:origin,
                            randomRopt ~ 0 + random.dummy:(1 | replicate),
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
                         seed = 3566+1,
                         cores = 4,
                         chains = 4,
                         file = "../models/growth_rate_model")

plot(growth_rate_model)
summary(growth_rate_model)

## Posterior predictive check
pdf("../figures/raw/fig_S6_ppc_growth_rate.pdf", height = 10, width = 10, pointsize = 3); set.seed(2491); pp_check(growth_rate_model, type = "ecdf_overlay_grouped", group = "replicate"); dev.off()

#### Plot growth rate TPCs ----

growth_rate_model_posteriors <- as.data.frame(growth_rate_model)

## TPCs with all data
pdf("../figures/raw/fig_S1_growth_rate_TPCs_all_data.pdf", height = 2, width = 2, pointsize = 3); set.seed(5285); for(i in c("anc", "cold", "hot")) {

  plot(NA,
       ylab = "Growth rate", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 0.00015), xlim = c(10, 50),
       main = i, cex.main = 1.75)

  for(g in unique(subset(dat, selection.regime == i)$origin)){

    tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

    points(jitter(growth.rate) ~ jitter(temperature),
           data = tempdat,
           pch = unique(tempdat$point), cex = 1,
           lwd = 0.5,
           bg = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[4],
           col = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[2])

    if(unique(tempdat$origin) == "bra"){linetype <- 1}
    if(unique(tempdat$origin) == "yem"){linetype <- 2}
    if(unique(tempdat$origin) == "ca"){linetype <- 3}

    TPC <- parse(text = paste0("curve(LRF(temp = x,
    Tmin = point_estimate(growth_rate_model_posteriors$`b_Tmin_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    Tmax = point_estimate(growth_rate_model_posteriors$`b_Tmax_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    Topt = point_estimate(growth_rate_model_posteriors$`b_Topt_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
      Ropt = point_estimate(growth_rate_model_posteriors$`b_Ropt_origin", g, ":selection.regime", i, "`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = linetype)"))

    eval(TPC)
  }

  groupmeans_regime <- aggregate(growth.rate ~ replicate + temperature + origin + colour + point, data = subset(dat, selection.regime == i & both.sexes == 1 & adult.offspring != 0), "geomean")
  groupmeans_regime <- aggregate(growth.rate ~ temperature + origin + colour + point, data = groupmeans_regime, "geomean")

  for(g in unique(groupmeans_regime$origin)){

    tempdat <- subset(groupmeans_regime, origin == g)

    points(growth.rate ~ temperature,
           data = tempdat,
           pch = unique(tempdat$point), cex = 2,
           lwd = 0.5,
           lty = 3,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  ytick <- c(0, 0.00005, 0.0001, 0.00015)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick * 1000, pos = 2, xpd = TRUE, cex=1.5)

}; legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n"); dev.off()

## Clean TPCs (without raw data)
pdf("../figures/raw/fig_S3_growth_rate_TPCs.pdf", height = 2, width = 2, pointsize = 3); set.seed(5285); for(i in c("anc", "cold", "hot")) {

  plot(NA,
       ylab = "Growth rate", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 0.00015), xlim = c(10, 50),
       main = i, cex.main = 1.75)

  abline(h = c(0, 0.00005, 0.0001, 0.00015), col = "grey40", lwd = 0.5, lty = 3)

  for(g in unique(subset(dat, selection.regime == i)$origin)){

    tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

    if(unique(tempdat$origin) == "bra"){linetype <- 1}
    if(unique(tempdat$origin) == "yem"){linetype <- 2}
    if(unique(tempdat$origin) == "ca"){linetype <- 3}

    TPC <- parse(text = paste0("curve(LRF(temp = x,
    Tmin = point_estimate(growth_rate_model_posteriors$`b_Tmin_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    Tmax = point_estimate(growth_rate_model_posteriors$`b_Tmax_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    Topt = point_estimate(growth_rate_model_posteriors$`b_Topt_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
      Ropt = point_estimate(growth_rate_model_posteriors$`b_Ropt_origin", g, ":selection.regime", i, "`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = linetype)"))

    eval(TPC)
  }

  groupmeans_regime <- aggregate(growth.rate ~ replicate + temperature + origin + colour + point, data = subset(dat, selection.regime == i & both.sexes == 1 & adult.offspring != 0), "geomean")
  groupmeans_regime <- aggregate(growth.rate ~ temperature + origin + colour + point, data = groupmeans_regime, "geomean")

  for(g in unique(groupmeans_regime$origin)){

    tempdat <- subset(groupmeans_regime, origin == g)

    points(growth.rate ~ temperature,
           data = tempdat,
           pch = unique(tempdat$point), cex = 2,
           lwd = 0.5,
           lty = 3,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  ytick <- c(0, 0.00005, 0.0001, 0.00015)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick * 1000, pos = 2, xpd = TRUE, cex=1.5)

}; legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n"); dev.off()

#### Plot posterior estimates ----

pdf("../figures/raw/fig_S2_growth_rate_posteriors.pdf", height = 2.5, width = 1.5, pointsize = 3); for(parameter in c("Tmin", "Topt", "Tmax", "Ropt")) {
  iter <- 1

  # Compute x-axis limits from 95% HDIs
  minima <- c()
  maxima <- c()
  for(i in c("anc", "cold", "hot")){
    for(g in unique(subset(dat, selection.regime == i)$origin)){
      eval(parse(text = paste0("minima <- c(minima, as.numeric(ci(growth_rate_model_posteriors$`b_", parameter, "_origin", g, ":selection.regime", i, "`, type ='HDI', 0.95)[2]))
                                                                    maxima <- c(maxima, as.numeric(ci(growth_rate_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, type ='HDI', 0.95)[3]))")))
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
        dens <- as.data.frame(density(growth_rate_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, n = 10000))
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

      MAP <- eval(parse(text = paste0("point_estimate(growth_rate_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, centrality ='MAP')")))
      hdpi <- eval(parse(text = paste0("ci(growth_rate_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, type ='HDI', 0.85)")))

      arrows(y0 = 50-ypos, y1 = 50-ypos,
             x0 = as.numeric(hdpi[2]), x1 = as.numeric(hdpi[3]),
             lwd = 2, code = 0,
             col = colorRampPalette(c(unique(subset(dat, selection.regime ==i)$colour), "black"))(6)[4])

      points(MAP, y = 50-ypos,
             pch =   unique(subset(dat, origin == g)$point),
             cex = 1.5,
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

#### Extract posterior estimates ----

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
trait <- "growth_rate"

## Main effects
for (parameter in parameters) {
  for (origin in unique(dat$origin)) {
    for (selection_regime in unique(dat$selection.regime)) {

      effect <- "main"
      param_name <- paste0("b_", parameter, "_origin", origin, ":selection.regime", selection_regime)

      if (param_name %in% names(growth_rate_model_posteriors)) {
        samples <- growth_rate_model_posteriors[[param_name]]
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

      if (fixed_param_name %in% names(growth_rate_model_posteriors) && random_param_name %in% names(growth_rate_model_posteriors)) {

        fixed_samples <- growth_rate_model_posteriors[[fixed_param_name]]
        random_samples <- growth_rate_model_posteriors[[random_param_name]]
        total_samples <- fixed_samples + random_samples

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

## Random effects (Ropt via randomRopt)
parameter <- "randomRopt"
for (replicate in unique(dat$replicate)) {

  random_dummy <- unique(dat$random.dummy[dat$replicate == replicate])

  if (random_dummy == 1) {

    effect <- "random"

    origin <- unique(dat$origin[dat$replicate == replicate])
    selection_regime <- unique(dat$selection.regime[dat$replicate == replicate])

    fixed_param_name <- paste0("b_Ropt_origin", origin, ":selection.regime", selection_regime)
    random_param_name <- paste0("r_replicate__", parameter, "[", replicate, ",Intercept]")

    if (fixed_param_name %in% names(growth_rate_model_posteriors) && random_param_name %in% names(growth_rate_model_posteriors)) {

      fixed_samples <- growth_rate_model_posteriors[[fixed_param_name]]
      random_samples <- growth_rate_model_posteriors[[random_param_name]]
      total_samples <- fixed_samples * exp(random_samples) # Ropt * exp(randomRopt)

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

posterior_estimates_growth_rate <- posterior_estimates
write.table(posterior_estimates_growth_rate, file = "../output/parameter_estimates_growth_rate.txt", sep = "\t", row.names = FALSE, quote = F)
