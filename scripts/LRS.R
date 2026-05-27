#### Load packages, functions, and data ----
source("packages_functions_data.R")

#### Oviposition rate / LRS ratio vs temperature ----

pdf("../figures/raw/fig_S16_oviposition_rate_LRS.pdf", height = 2.5, width = 2.5, pointsize = 3); {

  ## Raw data
  plot(NA, ylim = c(-4.9, -2.8), xlim = c(20, 38), ylab =" Log-transformed hourly\noviposition rate/LRS ratio", xlab = "Temperature (°C)", xaxt = "n")
  axis(1, at = c(23, 29, 35))
  for(i in unique(oviposition_rate_LRS$pop)){
  abline(h = c(coef(cor_fm2)[1], coef(cor_fm1)[1]), col = "orange", lwd = 1, lty = 3)
    if(unique(oviposition_rate_LRS$origin[which(oviposition_rate_LRS$pop == i)]) == "bra"){ char <- 21 }
    if(unique(oviposition_rate_LRS$origin[which(oviposition_rate_LRS$pop == i)]) == "yem"){ char <- 22 }
    if(unique(oviposition_rate_LRS$origin[which(oviposition_rate_LRS$pop == i)]) == "ca"){ char <- 24 }
    lines(log.ratio ~ temperature, lwd = 1, pch = 21,
          col = c("darkgrey", "dodgerblue", "orangered")[which(unique(oviposition_rate_LRS$selection.regime)==unique(oviposition_rate_LRS[oviposition_rate_LRS$pop == i,]$selection.regime))],
          data = oviposition_rate_LRS[oviposition_rate_LRS$pop == i,])
    points(log.ratio ~ temperature, cex = 1.5, lwd = 0.5, pch = char,
           bg = c("darkgrey", "dodgerblue", "orangered")[which(unique(oviposition_rate_LRS$selection.regime)==unique(oviposition_rate_LRS[oviposition_rate_LRS$pop == i,]$selection.regime))],
           col = "black", data = oviposition_rate_LRS[oviposition_rate_LRS$pop == i,])
  }

  ## Linear model
  fm <- lm(log.ratio ~ 0 + selection.regime:origin + I(temperature-23):selection.regime, data = oviposition_rate_LRS)
  summary(fm)

  coef_summary <- summary(fm)$coefficients
  point_estimates <- coef_summary[, "Estimate"]
  se <- coef_summary[, "Std. Error"]

  # 85% confidence intervals
  lower_ci <- point_estimates - 1.440 * se
  upper_ci <- point_estimates + 1.440 * se

  ## Intercept plot
  plot(NA, ylim = c(exp(-4.9), exp(-3.6)), xlim = c(0.5, 3.5), ylab = "Intercept\n(Hourly oviposition rate/LRS ratio at 23°C)", xlab = "Selection regime", xaxt = "n")
  axis(1, at = c(1, 2, 3), labels = c("Ancestral", "Cold", "Hot"))
  counter = 0
  for(i in 1:3){
    color <- c("darkgrey", "dodgerblue", "orangered")[i]
    sr <- c("selection.regimeanc:originbra", "selection.regimeanc:originyem", "selection.regimeanc:originca",
            "selection.regimec:originbra", "selection.regimec:originyem", "selection.regimec:originca",
            "selection.regimeh:originbra", "selection.regimeh:originyem", "selection.regimeh:originca")[i + counter]

    arrows(x0 = i - 0.25, x1 = i - 0.25,
           y0 = exp(lower_ci[which(names(lower_ci) == sr)]),
           y1 = exp(upper_ci[which(names(lower_ci) == sr)]), code = 3, angle = 90, length = 0.05, lwd = 1, col = color)
    points(x = i - 0.25, y = exp(point_estimates[which(names(lower_ci) == sr)]), bg = color, pch = 21, cex = 1.5, lwd = 0.5)

    sr <- c("selection.regimeanc:originbra", "selection.regimeanc:originyem", "selection.regimeanc:originca",
            "selection.regimec:originbra", "selection.regimec:originyem", "selection.regimec:originca",
            "selection.regimeh:originbra", "selection.regimeh:originyem", "selection.regimeh:originca")[i + counter + 1]

    arrows(x0 = i, x1 = i,
           y0 = exp(lower_ci[which(names(lower_ci) == sr)]),
           y1 = exp(upper_ci[which(names(lower_ci) == sr)]), code = 3, angle = 90, length = 0.05, lwd = 1, col = color)
    points(x = i, y = exp(point_estimates[which(names(lower_ci) == sr)]), bg = color, pch = 22, cex = 1.5, lwd = 0.5)

    sr <- c("selection.regimeanc:originbra", "selection.regimeanc:originyem", "selection.regimeanc:originca",
            "selection.regimec:originbra", "selection.regimec:originyem", "selection.regimec:originca",
            "selection.regimeh:originbra", "selection.regimeh:originyem", "selection.regimeh:originca")[i + counter + 2]

    arrows(x0 = i + 0.25, x1 = i + 0.25,
           y0 = exp(lower_ci[which(names(lower_ci) == sr)]),
           y1 = exp(upper_ci[which(names(lower_ci) == sr)]), code = 3, angle = 90, length = 0.05, lwd = 1, col = color)
    points(x = i + 0.25, y = exp(point_estimates[which(names(lower_ci) == sr)]), bg = color, pch = 24, cex = 1.5, lwd = 0.5)
    counter = counter + 2
  }

  ## Slope plot
  plot(NA, ylim = c(0.03, 0.08), xlim = c(0.5, 3.5), ylab = "Slope\n(Hourly oviposition rate/LRS ratio ~ Temperature)", xlab = "Selection regime", xaxt = "n")
  axis(1, at = c(1, 2, 3), labels = c("Ancestral", "Cold", "Hot"))
  for(i in 1:3){
    sr <- c("selection.regimeanc:I(temperature - 23)", "selection.regimec:I(temperature - 23)", "selection.regimeh:I(temperature - 23)")[i]
    color <- c("darkgrey", "dodgerblue", "orangered")[i]

    arrows(x0 = i, x1 = i,
           y0 = lower_ci[which(names(lower_ci) == sr)],
           y1 = upper_ci[which(names(lower_ci) == sr)], code = 3, angle = 90, length = 0.05, lwd = 1, col = color)
    points(x = i, y = point_estimates[which(names(lower_ci) == sr)], bg = color, pch = 21, cex = 1.5, lwd = 0.5)
  }
}; dev.off()

#### Oviposition rate / LRS ratio vs temperature (batch-corrected) ----

pdf("../figures/raw/fig_S16_oviposition_rate_LRS_corrected.pdf", height = 2.5, width = 2.5, pointsize = 3); {

  ## Raw data
  plot(NA, ylim = c(-4.9, -2.8), xlim = c(20, 38), ylab =" Log-transformed hourly\noviposition rate/LRS ratio", xlab = "Temperature (°C)", xaxt = "n")
  axis(1, at = c(23, 29, 35))
  for(i in unique(oviposition_rate_LRS$pop)){
    if(unique(oviposition_rate_LRS$origin[which(oviposition_rate_LRS$pop == i)]) == "bra"){ char <- 21 }
    if(unique(oviposition_rate_LRS$origin[which(oviposition_rate_LRS$pop == i)]) == "yem"){ char <- 22 }
    if(unique(oviposition_rate_LRS$origin[which(oviposition_rate_LRS$pop == i)]) == "ca"){ char <- 24 }
    lines(log.ratio.corrected ~ temperature, lwd = 1, pch = 21,
          col = c("darkgrey", "dodgerblue", "orangered")[which(unique(oviposition_rate_LRS$selection.regime)==unique(oviposition_rate_LRS[oviposition_rate_LRS$pop == i,]$selection.regime))],
          data = oviposition_rate_LRS[oviposition_rate_LRS$pop == i,])
    points(log.ratio.corrected ~ temperature, cex = 1.5, lwd = 0.5, pch = char,
           bg = c("darkgrey", "dodgerblue", "orangered")[which(unique(oviposition_rate_LRS$selection.regime)==unique(oviposition_rate_LRS[oviposition_rate_LRS$pop == i,]$selection.regime))],
           col = "black", data = oviposition_rate_LRS[oviposition_rate_LRS$pop == i,])
  }

  ## Linear model
  fm <- lm(log.ratio.corrected ~ 0 + selection.regime:origin + I(temperature-23):selection.regime, data = oviposition_rate_LRS)
  summary(fm)

  coef_summary <- summary(fm)$coefficients
  point_estimates <- coef_summary[, "Estimate"]
  se <- coef_summary[, "Std. Error"]

  # 85% confidence intervals
  lower_ci <- point_estimates - 1.440 * se
  upper_ci <- point_estimates + 1.440 * se

  ## Intercept plot
  plot(NA, ylim = c(exp(-4.9), exp(-3.6)), xlim = c(0.5, 3.5), ylab = "Intercept\n(Hourly oviposition rate/LRS ratio at 23°C)", xlab = "Selection regime", xaxt = "n")
  axis(1, at = c(1, 2, 3), labels = c("Ancestral", "Cold", "Hot"))
  counter = 0
  for(i in 1:3){
    color <- c("darkgrey", "dodgerblue", "orangered")[i]
    sr <- c("selection.regimeanc:originbra", "selection.regimeanc:originyem", "selection.regimeanc:originca",
            "selection.regimec:originbra", "selection.regimec:originyem", "selection.regimec:originca",
            "selection.regimeh:originbra", "selection.regimeh:originyem", "selection.regimeh:originca")[i + counter]

    arrows(x0 = i - 0.25, x1 = i - 0.25,
           y0 = exp(lower_ci[which(names(lower_ci) == sr)]),
           y1 = exp(upper_ci[which(names(lower_ci) == sr)]), code = 3, angle = 90, length = 0.05, lwd = 1, col = color)
    points(x = i - 0.25, y = exp(point_estimates[which(names(lower_ci) == sr)]), bg = color, pch = 21, cex = 1.5, lwd = 0.5)

    sr <- c("selection.regimeanc:originbra", "selection.regimeanc:originyem", "selection.regimeanc:originca",
            "selection.regimec:originbra", "selection.regimec:originyem", "selection.regimec:originca",
            "selection.regimeh:originbra", "selection.regimeh:originyem", "selection.regimeh:originca")[i + counter + 1]

    arrows(x0 = i, x1 = i,
           y0 = exp(lower_ci[which(names(lower_ci) == sr)]),
           y1 = exp(upper_ci[which(names(lower_ci) == sr)]), code = 3, angle = 90, length = 0.05, lwd = 1, col = color)
    points(x = i, y = exp(point_estimates[which(names(lower_ci) == sr)]), bg = color, pch = 22, cex = 1.5, lwd = 0.5)

    sr <- c("selection.regimeanc:originbra", "selection.regimeanc:originyem", "selection.regimeanc:originca",
            "selection.regimec:originbra", "selection.regimec:originyem", "selection.regimec:originca",
            "selection.regimeh:originbra", "selection.regimeh:originyem", "selection.regimeh:originca")[i + counter + 2]

    arrows(x0 = i + 0.25, x1 = i + 0.25,
           y0 = exp(lower_ci[which(names(lower_ci) == sr)]),
           y1 = exp(upper_ci[which(names(lower_ci) == sr)]), code = 3, angle = 90, length = 0.05, lwd = 1, col = color)
    points(x = i + 0.25, y = exp(point_estimates[which(names(lower_ci) == sr)]), bg = color, pch = 24, cex = 1.5, lwd = 0.5)
    counter = counter + 2
  }

  ## Slope plot
  plot(NA, ylim = c(0.03, 0.08), xlim = c(0.5, 3.5), ylab = "Slope\n(Hourly oviposition rate/LRS ratio ~ Temperature)", xlab = "Selection regime", xaxt = "n")
  axis(1, at = c(1, 2, 3), labels = c("Ancestral", "Cold", "Hot"))
  for(i in 1:3){
    sr <- c("selection.regimeanc:I(temperature - 23)", "selection.regimec:I(temperature - 23)", "selection.regimeh:I(temperature - 23)")[i]
    color <- c("darkgrey", "dodgerblue", "orangered")[i]

    arrows(x0 = i, x1 = i,
           y0 = lower_ci[which(names(lower_ci) == sr)],
           y1 = upper_ci[which(names(lower_ci) == sr)], code = 3, angle = 90, length = 0.05, lwd = 1, col = color)
    points(x = i, y = point_estimates[which(names(lower_ci) == sr)], bg = color, pch = 21, cex = 1.5, lwd = 0.5)
  }
}; dev.off()

#### Calculate oviposition rates ----

dat$oviposition.rate <- ifelse(dat$selection.regime == "anc",
                               oviposition_rate_function_anc_hourly(temp = dat$temperature, LRS = dat$adult.offspring),
                               oviposition_rate_function_selection_hourly(temp = dat$temperature, LRS = dat$adult.offspring))

#### Plot raw fecundity data ----

pdf("../figures/raw/fig_S17_fecundity_data.pdf", height = 2, width = 2, pointsize = 3); for(i in c("anc", "cold", "hot")) {

  plot(NA,
       ylab = "Fecundity", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 120), xlim = c(10, 50),
       main = i, cex.main = 1.75)
  abline(h = c(0, 60, 120), lty = 3, col = "grey")

  for(g in unique(subset(dat, selection.regime == i)$origin)){
    tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)
    points(jitter(adult.offspring) ~ jitter(temperature),
           data = tempdat,
           pch = unique(tempdat$point), cex = 1,
           lwd = 0.5,
           bg = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[4],
           col = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[2])
  }

  ## Group means
  groupmeans_regime <- aggregate(adult.offspring ~  temperature + origin + colour + point, data = subset(dat, selection.regime == i & both.sexes == 1 & adult.offspring != 0), "mean")

  for(g in unique(groupmeans_regime$origin)){
    tempdat <- subset(groupmeans_regime, origin == g)
    if(unique(tempdat$origin) == "bra"){linetype <- 1}
    if(unique(tempdat$origin) == "yem"){linetype <- 2}
    if(unique(tempdat$origin) == "ca"){linetype <- 3}
    lines(adult.offspring ~ temperature,
          data = tempdat,
          lwd = 0.5,
          lty = linetype,
          col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  for(g in unique(groupmeans_regime$origin)){
    tempdat <- subset(groupmeans_regime, origin == g)
    if(unique(tempdat$origin) == "bra"){point <- 21}
    if(unique(tempdat$origin) == "yem"){point <- 22}
    if(unique(tempdat$origin) == "ca"){point <- 24}
    points(adult.offspring ~ temperature,
           data = tempdat,
           pch = point, cex = 2,
           lwd = 0.5,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  ytick <- c(0, 60, 120)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)

}; legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n"); dev.off()

#### Plot oviposition rate data ----

pdf("../figures/raw/fig_S17_oviposition_rate_data.pdf", height = 2, width = 2, pointsize = 3); for(i in c("anc", "cold","hot")) {

  plot(NA,
       ylab = "Oviposition rate", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 2), xlim = c(10, 50),
       main = i, cex.main = 1.75)
  abline(h = c(0, 1, 2), lty = 3, col = "grey")

  for(g in unique(subset(dat, selection.regime == i)$origin)){
    tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)
    points(jitter(oviposition.rate) ~ jitter(temperature),
           data = tempdat,
           pch = unique(tempdat$point), cex = 1,
           lwd = 0.5,
           bg = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[4],
           col = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[2])
  }

  ## Group means
  groupmeans_regime <- aggregate(oviposition.rate ~  temperature + origin + colour + point, data = subset(dat, selection.regime == i & both.sexes == 1 & oviposition.rate != 0), "mean")

  for(g in unique(groupmeans_regime$origin)){
    tempdat <- subset(groupmeans_regime, origin == g)
    if(unique(tempdat$origin) == "bra"){linetype <- 1}
    if(unique(tempdat$origin) == "yem"){linetype <- 2}
    if(unique(tempdat$origin) == "ca"){linetype <- 3}
    lines(oviposition.rate ~ temperature,
          data = tempdat,
          lwd = 0.5,
          lty = linetype,
          col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  for(g in unique(groupmeans_regime$origin)){
    tempdat <- subset(groupmeans_regime, origin == g)
    if(unique(tempdat$origin) == "bra"){point <- 21}
    if(unique(tempdat$origin) == "yem"){point <- 22}
    if(unique(tempdat$origin) == "ca"){point <- 24}
    points(oviposition.rate ~ temperature,
           data = tempdat,
           pch = point, cex = 2,
           lwd = 0.5,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  ytick <- c(0, 1, 2)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)

}; dev.off()

#### Bayesian fecundity model ----

## Priors
fecundity_model_prior <- prior(normal(29, 5), nlpar="Topt", lb = 0) +
  prior(normal(70,25), nlpar="fecmax", lb = 0) +
  prior(lognormal(log(22), log(1.1)), nlpar="breadth", lb = 0) +
  prior(lognormal((2.5),(2.5)), nlpar="fatness", lb = 0) +
  prior(lognormal(log(0.03), log(1.2)), nlpar="a", lb = 0) +
  prior(normal(28, 2), nlpar="h") +
  prior(normal(-4.5, 2), nlpar="k")+
  prior(normal(0, 2.5), dpar="shape")

## Starting values
init <- parse(text=paste0("list(b_Topt = as.array(rnorm(9, 29, 1.5)),
  b_fecmax = as.array(rnorm(9, 70, 1.5)),
  b_breadth = as.array(rlnorm(9, log(25), log(1.1))),
  b_fatness = as.array(rlnorm(2, log(4), log(1.5))))"))
set.seed(1646);inits_list <- list(eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init))

## Fit model
offspring_model <- brm(bf(adult.offspring ~ log(fec_function(temperature, Topt, fecmax, breadth, fatness)),
                          nlf(zi ~ a * (temperature - h)^2 + k),
                          shape ~ 0 + temperature.factor:selection.regime,
                          fatness  ~ 0 + anc.dummy,
                          a + h + k ~ 0 + selection.regime,
                          Topt + fecmax + breadth  ~ 0 + origin:selection.regime + random.dummy:(1 | replicate),
                          nl = T),
                       data = dat,
                       family = zero_inflated_negbinomial,
                       prior = fecundity_model_prior,
                       control = list(adapt_delta = 0.9, max_treedepth = 15),
                       warmup = 2000,
                       init = inits_list,
                       stanvars = stanvars,
                       iter = 3000,
                       thin = 1,
                       seed = 5559,
                       cores = 6,
                       chains = 6,
                       file = "../models/offspring_model")

# plot(offspring_model)
summary(offspring_model)

## Posterior predictive check
pdf("../figures/raw/fig_S5_ppc_LRS.pdf", height = 10, width = 10, pointsize = 3); set.seed(2491); print(pp_check(offspring_model, type = "ecdf_overlay_grouped", group = "replicate")); dev.off()

## Extract posteriors
offspring_model_posteriors <- as.data.frame(offspring_model)
colnames(offspring_model_posteriors)

#### Plot fecundity and non-viability curves ----

pdf("../figures/raw/fig_S1_fecundity_curves.pdf", height = 2, width = 2, pointsize = 3); set.seed(1589); for(i in c("anc", "cold", "hot")) {

  plot(NA,
       ylab = "Fecundity", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 120), xlim = c(10, 50),
       main = i, cex.main = 1.75)

  for(g in unique(subset(dat, selection.regime == i)$origin)){
    tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

    points(jitter(adult.offspring) ~ jitter(temperature),
           data = tempdat,
           pch = unique(tempdat$point), cex = 1,
           lwd = 0.5,
           bg = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[4],
           col = colorRampPalette(c(unique(tempdat$colour), "white"))(6)[2])

    if(unique(tempdat$origin) == "bra"){linetype <- 1}
    if(unique(tempdat$origin) == "yem"){linetype <- 2}
    if(unique(tempdat$origin) == "ca"){linetype <- 3}
    if(i == "anc"){
    fec <- parse(text = paste0("curve(fec_function(temp = x,
    Topt = point_estimate(offspring_model_posteriors$`b_Topt_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    fecmax = point_estimate(offspring_model_posteriors$`b_fecmax_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    breadth = point_estimate(offspring_model_posteriors$`b_breadth_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
      fatness = point_estimate(offspring_model_posteriors$`b_fatness_anc.dummyanc`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = linetype)"))
    }
    else{
      fec <- parse(text = paste0("curve(fec_function(temp = x,
    Topt = point_estimate(offspring_model_posteriors$`b_Topt_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    fecmax = point_estimate(offspring_model_posteriors$`b_fecmax_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    breadth = point_estimate(offspring_model_posteriors$`b_breadth_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
      fatness = point_estimate(offspring_model_posteriors$`b_fatness_anc.dummynot`, centrality ='MAP')),
      10, 50, add = T, col = unique(tempdat$colour), lty = linetype)"))
    }
    eval(fec)

    viability <- parse(text = paste0("curve(invlogit(viab_function(temp = x,
    a = point_estimate(offspring_model_posteriors$`b_a_selection.regime", i, "`, centrality ='MAP'),
    h = point_estimate(offspring_model_posteriors$`b_h_selection.regime", i, "`, centrality ='MAP'),
    k = point_estimate(offspring_model_posteriors$`b_k_selection.regime", i, "`, centrality ='MAP'))) * 120,
      10, 50, add = T, col = colorRampPalette(c(unique(tempdat$colour), 'black'))(6)[4], lty = 1, lwd = 1)"))
    eval(viability)
  }

  ## Group means
  groupmeans_regime <- aggregate(adult.offspring ~  temperature + origin + colour + point, data = subset(dat, selection.regime == i & both.sexes == 1 & adult.offspring != 0), "mean")

  for(g in unique(groupmeans_regime$origin)){
    tempdat <- subset(groupmeans_regime, origin == g)
    if(unique(tempdat$origin) == "bra"){point <- 21}
    if(unique(tempdat$origin) == "yem"){point <- 22}
    if(unique(tempdat$origin) == "ca"){point <- 24}
    points(adult.offspring ~ temperature,
           data = tempdat,
           pch = point, cex = 2,
           lwd = 0.5,
           lty = 3,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  ytick <- c(0, 60, 120)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)

  }; legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n"); dev.off()

#### Plot oviposition rate curves ----

pdf("../figures/raw/fig_S3_oviposition_rate_curves.pdf", height = 2, width = 2, pointsize = 3); for(i in c("anc", "cold", "hot")) {

  plot(NA,
       ylab = "Oviposition rate (eggs / hour)", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 40), xlim = c(10, 50),
       main = i, cex.main = 1.75)
  abline(h = c(0, 20, 40), col = "grey40", lwd = 0.5, lty = 3)

  for(g in unique(subset(dat, selection.regime == i)$origin)){
    tempdat <- subset(dat, selection.regime == i & both.sexes == 1 & origin == g)

    if(unique(tempdat$origin) == "bra"){linetype <- 1}
    if(unique(tempdat$origin) == "yem"){linetype <- 2}
    if(unique(tempdat$origin) == "ca"){linetype <- 3}
    if(i == "anc"){
    fec <- parse(text = paste0("curve(oviposition_rate_function_anc(temp = x, LRS = fec_function(temp = x,
    Topt = point_estimate(offspring_model_posteriors$`b_Topt_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    fecmax = point_estimate(offspring_model_posteriors$`b_fecmax_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    breadth = point_estimate(offspring_model_posteriors$`b_breadth_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
      fatness = point_estimate(offspring_model_posteriors$`b_fatness_anc.dummyanc`, centrality ='MAP'))),
      10, 50, add = T, col = unique(tempdat$colour), lty = linetype)"))
    }
    else{
      fec <- parse(text = paste0("curve(oviposition_rate_function_selection(temp = x, LRS = fec_function(temp = x,
    Topt = point_estimate(offspring_model_posteriors$`b_Topt_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    fecmax = point_estimate(offspring_model_posteriors$`b_fecmax_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
    breadth = point_estimate(offspring_model_posteriors$`b_breadth_origin", g, ":selection.regime", i, "`, centrality ='MAP'),
      fatness = point_estimate(offspring_model_posteriors$`b_fatness_anc.dummynot`, centrality ='MAP'))),
      10, 50, add = T, col = unique(tempdat$colour), lty = linetype)"))
    }
    eval(fec)
  }

  ## Group means
  groupmeans_regime <- aggregate(oviposition.rate ~  temperature + origin + colour + point, data = subset(dat, selection.regime == i & both.sexes == 1 & adult.offspring != 0), "mean")

  for(g in unique(groupmeans_regime$origin)){
    tempdat <- subset(groupmeans_regime, origin == g)
    if(unique(tempdat$origin) == "bra"){point <- 21}
    if(unique(tempdat$origin) == "yem"){point <- 22}
    if(unique(tempdat$origin) == "ca"){point <- 24}
    points(oviposition.rate * 24 ~ temperature,
           data = tempdat,
           pch = point, cex = 2,
           lwd = 0.5,
           lty = 3,
           bg = unique(tempdat$colour),
           col = colorRampPalette(c(unique(tempdat$colour), "black"))(6)[4])
  }

  ytick <- c(0, 20, 40)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick,  par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)

}; legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n"); dev.off()

#### Plot viability rate curves ----

parameter_estimates_dev_rate  <- read.delim("../output/parameter_estimates_dev_rate.txt")
parameter_estimates_viab_fec  <- read.delim("../output/parameter_estimates_viab_fec.txt")

pdf("../figures/raw/fig_S3_viability_curves.pdf", height = 2, width = 2, pointsize = 3); for(i in c("anc", "cold", "hot")) {

  plot(NA,
       ylab = "Viability rate", xlab = "Temperature (°C)", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(0, 1), xlim = c(10, 50),
       main = i, cex.main = 1.75)
  abline(h = seq(0, 1, by = 0.25), col = "grey40", lwd = 0.5, lty = 3)

  for(g in unique(subset(dat, selection.regime == i)$origin)){
    tempdat <- subset(dat, selection.regime == i & origin == g)

    if(unique(tempdat$origin) == "bra"){linetype <- 1; point <- 21}
    if(unique(tempdat$origin) == "yem"){linetype <- 2; point <- 22}
    if(unique(tempdat$origin) == "ca"){linetype <- 3; point <- 24}

    colour_map <- c("anc" = "darkgrey", "cold" = "dodgerblue", "hot" = "orangered")
    colour <- colour_map[i]

    # Development rate parameters
    dev_params <- parameter_estimates_dev_rate %>%
      filter(trait == "development_rate",
             effect == "main",
             origin == g,
             selection.regime == i,
             parameter %in% c("Tmin", "Tmax", "Topt", "Ropt")) %>%
      select(parameter, MAP) %>%
      pivot_wider(names_from = parameter, values_from = MAP)

    Tmin <- as.numeric(dev_params$Tmin)
    Tmax <- as.numeric(dev_params$Tmax)
    Topt <- as.numeric(dev_params$Topt)
    Ropt <- as.numeric(dev_params$Ropt)

    # Viability parameters (per selection regime, not per origin)
    viab_params <- parameter_estimates_viab_fec %>%
      filter(trait == "viability",
             effect == "main",
             selection.regime == i,
             parameter %in% c("a", "h", "k")) %>%
      select(parameter, MAP) %>%
      pivot_wider(names_from = parameter, values_from = MAP)

    a <- as.numeric(viab_params$a)
    h <- as.numeric(viab_params$h)
    k <- as.numeric(viab_params$k)

    curve(viability_rate(x, Tmin, Tmax, Topt, Ropt, a, h, k),
          from = 10, to = 50, col = colour, lty = linetype, add = TRUE)
  }

  ytick <- seq(0, 1, by = 0.25)
  axis(side=2, at=ytick, labels = FALSE)
  text(y=ytick, par("usr")[1],
       labels = ytick, pos = 2, xpd = TRUE, cex=1.5)
}; legend("topleft", c("Brazil", "Yemen", "California"), pch = c(21, 22, 24), pt.bg = "grey", pt.lwd = 0.5, lty = c(1, 2, 3), cex = 1.5, bty = "n"); dev.off()

#### Plot posterior estimates ----

pdf("../figures/raw/fig_S2_fecundity_posteriors.pdf", height = 2.5, width = 1.5, pointsize = 3); for(parameter in c("Topt", "fecmax", "breadth", "h", "k", "a")){

    if(parameter == "Topt" | parameter == "fecmax" | parameter == "breadth"){

      minima <- c()
      maxima <- c()
      for(i in c("anc", "cold", "hot")){
        for(g in unique(subset(dat, selection.regime == i)$origin)){
          eval(parse(text = paste0("minima <- c(minima, as.numeric(ci(offspring_model_posteriors$`b_", parameter, "_origin", g, ":selection.regime", i, "`, type ='HDI', 0.95)[2]))
                                                                    maxima <- c(maxima, as.numeric(ci(offspring_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, type ='HDI', 0.95)[3]))")))
        }
      }

  ## Per-origin posteriors
  iter <- 1
  ypos <- iter
  mean_estimate <- c()

  plot(NA,
       ylab = "", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
       ylim = c(37, 51), xlim = c(min(minima), max(maxima)),
       main = parameter, cex.main = 1.75)

  for(i in c("anc", "cold", "hot")){
    for(g in unique(subset(dat, selection.regime == i)$origin)){

      eval(parse(text = paste0("{
        dens <- as.data.frame(density(offspring_model_posteriors$`b_", parameter,"_origin", g, ":selection.regime", i, "`, n = 10000))
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

      MAP <- eval(parse(text = paste0("point_estimate(na.omit(offspring_model_posteriors$`b_", parameter, "_origin", g, ":selection.regime", i, "`), centrality ='MAP')")))
      hdpi <- eval(parse(text = paste0("ci(na.omit(offspring_model_posteriors$`b_", parameter, "_origin", g, ":selection.regime", i, "`), type ='HDI', 0.85)")))

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

  }

  else{

    minima <- c()
    maxima <- c()
    for(i in c("anc", "cold", "hot")){
      for(g in unique(subset(dat, selection.regime == i)$origin)){
        eval(parse(text = paste0("minima <- c(minima, as.numeric(ci(offspring_model_posteriors$`b_", parameter, "_selection.regime", i, "`, type ='HDI', 0.95)[2]))
                                                                    maxima <- c(maxima, as.numeric(ci(offspring_model_posteriors$`b_", parameter, "_selection.regime", i, "`, type ='HDI', 0.95)[3]))")))
      }
    }

    ## Per-regime posteriors
    iter <- 1
    ypos <- iter+1
    mean_estimate <- c()

    plot(NA,
         ylab = "", xlab = "Temperature", cex.axis = 1.5, cex.lab = 1.5, yaxt = "n",
         ylim = c(37, 51), xlim = c(min(minima), max(maxima)),
         main = parameter, cex.main = 1.75)

    for(i in c("anc", "cold", "hot")){

      eval(parse(text = paste0("{
        dens <- as.data.frame(density(offspring_model_posteriors$`b_", parameter, "_selection.regime", i, "`, n = 10000))
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

      MAP <- eval(parse(text = paste0("point_estimate(offspring_model_posteriors$`b_",  parameter, "_selection.regime", i, "`, centrality ='MAP')")))
      hdpi <- eval(parse(text = paste0("ci(offspring_model_posteriors$`b_", parameter,"_selection.regime", i, "`, type ='HDI', 0.85)")))

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

      ypos = ypos + 3
      mean_estimate <- c(mean_estimate, as.numeric(MAP))

      ypos = ypos + 1.5
      iter = iter + 1
    }

    abline(v = mean(mean_estimate), lty = 3)

  }
  abline(v = mean(mean_estimate), lty = 3)

}; dev.off()

#### Extract and save posterior estimates ----

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
trait_fecundity <- "fecundity"

for (parameter in fecundity_parameters) {
  for (origin in unique(dat$origin)) {
    for (selection_regime in unique(dat$selection.regime)) {

      effect <- "main"
      param_name <- paste0("b_", parameter, "_origin", origin, ":selection.regime", selection_regime)

      if (parameter == "fatness"){
        if(selection_regime == "anc"){
        samples <- offspring_model_posteriors[["b_fatness_anc.dummyanc"]]
        }
        else{
          samples <- offspring_model_posteriors[["b_fatness_anc.dummynot"]]
        }
        MAP <- point_estimate(samples, centrality='MAP')
        HPDI <- ci(samples, method='HDI', ci=0.9)

        posterior_estimates <- rbind(posterior_estimates,
                                     data.frame(trait=trait_fecundity,
                                                parameter=parameter,
                                                effect=effect,
                                                origin=origin,
                                                selection.regime=selection_regime,
                                                replicate=NA,
                                                MAP=MAP,
                                                hpdi.lower=HPDI$CI_low,
                                                hpdi.upper=HPDI$CI_high,
                                                stringsAsFactors=FALSE))

      }

      else if (param_name %in% names(offspring_model_posteriors)) {
        samples <- offspring_model_posteriors[[param_name]]
        MAP <- point_estimate(samples, centrality='MAP')
        HPDI <- ci(samples, method='HDI', ci=0.9)

        posterior_estimates <- rbind(posterior_estimates,
                                     data.frame(trait=trait_fecundity,
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

## Random effects for replicates
for (parameter in fecundity_parameters[-4]) {
  for (replicate in unique(dat$replicate)) {

    random_dummy <- unique(dat$random.dummy[dat$replicate == replicate])

    if (random_dummy == 1) {

      effect <- "random"
      origin <- unique(dat$origin[dat$replicate == replicate])
      selection_regime <- unique(dat$selection.regime[dat$replicate == replicate])

      fixed_param_name <- paste0("b_", parameter, "_origin", origin, ":selection.regime", selection_regime)
      random_param_name <- paste0("r_replicate__", parameter, "[", replicate, ",Intercept]")

      if (fixed_param_name %in% names(offspring_model_posteriors) && random_param_name %in% names(offspring_model_posteriors)) {

        fixed_samples <- offspring_model_posteriors[[fixed_param_name]]
        random_samples <- offspring_model_posteriors[[random_param_name]]
        total_samples <- fixed_samples + random_samples

        MAP <- point_estimate(total_samples, centrality='MAP')
        HPDI <- ci(total_samples, method='HDI', ci=0.9)

        posterior_estimates <- rbind(posterior_estimates,
                                     data.frame(trait=trait_fecundity,
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

## Viability parameters
viability_parameters <- c("h", "k", "a")
trait_viability <- "viability"

for (parameter in viability_parameters) {
  for (selection_regime in unique(dat$selection.regime)) {
    effect <- "main"
    param_name <- paste0("b_", parameter, "_selection.regime", selection_regime)

    if (param_name %in% names(offspring_model_posteriors)) {
      samples <- offspring_model_posteriors[[param_name]]
      MAP <- point_estimate(samples, centrality='MAP')
      HPDI <- ci(samples, method='HDI', ci=0.9)

      posterior_estimates <- rbind(posterior_estimates,
                                   data.frame(trait=trait_viability,
                                              parameter=parameter,
                                              effect=effect,
                                              origin=NA,
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

## Save estimates
posterior_estimates_viab_fec <- posterior_estimates
write.table(posterior_estimates_viab_fec, file = "../output/parameter_estimates_viab_fec.txt", sep = "\t", row.names = FALSE, quote = F)
