# function to simulate the data set for fitting linear mixed model
simudata.lmm <- function (params) {
  # This function simulates a dataset for fitting linear mixed mode.
  # A more detailed explanation could be found [here](https://haiyangjin.github.io/portfolio/simulate-data-lmm/)
  # author: Haiyang Jin (https://haiyangjin.github.io/)

  # # HOW to use this function?
  # source("path/to/simulate_data_lmm.R")
  # params <- list(
  #   # set the parameters for experiment design
  #   IV1.levels = c("upright", "inverted"),  # (has to be two levels)
  #   IV2.levels = c("intact", "scrambled"),  # (has to be two levels)
  #   num_Subj = 30,
  #   num_Stim = 20,
  # 
  #   # set the mu for every bin (every condition)
  #   IV1.1_IV2.1 = 500,
  #   IV1.2_IV2.1 = 600,
  #   IV1.1_IV2.2 = 800,
  #   IV1.2_IV2.2 = 850,
  #
  #   # set the variances for lmm (std)
  #   var_residual = 30,  # residual
  #   var_rnd_int_subj = 40,  # random intercept for Subject
  #   var_rnd_int_stim = 50,  # random intercept for Stimuli
  #   var_rnd_slp_IV1_subj = 60,  # random slope of IV1 for Subject
  #   var_rnd_slp_IV2_subj = 70,  # random slope of IV2 for Subject
  #   var_rnd_slp_inter_subj = 20,  # random slope of IV1*IV2 for Subject
  #   var_rnd_slp_IV1_stim = 80,  # random slope of IV1 for Stimuli
  #   var_rnd_slp_IV2_stim = 90,  # random slope of IV2 for Stimuli
  #   var_rnd_slp_inter_stim = 20  # random slope of IV1*IV2 for Stimuli
  # )
  #
  # Usage: simudata <- simudata.lmm(params = params)
  
  # load library
  

  ## load the parameters from params
  # set the parameters for experiment design
  IV1.levels <- params$IV1.levels
  IV2.levels <- params$IV2.levels
  num_Subj <- params$num_Subj
  num_Stim <- params$num_Stim

  # set the mu for every bin (every condition)
  IV1.1_IV2.1 <- params$IV1.1_IV2.1
  IV1.2_IV2.1 <- params$IV1.2_IV2.1
  IV1.1_IV2.2 <- params$IV1.1_IV2.2
  IV1.2_IV2.2 <- params$IV1.2_IV2.2

  # set the variances for lmm (std)
  var_residual <- params$var_residual  # residual
  var_rnd_int_subj <- params$var_rnd_int_subj  # random intercept for Subject
  var_rnd_int_stim <- params$var_rnd_int_stim  # random intercept for Stimuli
  var_rnd_slp_IV1_subj <- params$var_rnd_slp_IV1_subj  # random slope of IV1 for Subject
  var_rnd_slp_IV2_subj <- params$var_rnd_slp_IV2_subj  # random slope of IV2 for Subject
  var_rnd_slp_inter_subj <- params$var_rnd_slp_inter_subj  # random slope of IV1*IV2 for Subject
  var_rnd_slp_IV1_stim <- params$var_rnd_slp_IV1_stim  # random slope of IV1 for Stimuli
  var_rnd_slp_IV2_stim <- params$var_rnd_slp_IV2_stim  # random slope of IV2 for Stimuli
  var_rnd_slp_inter_stim <- params$var_rnd_slp_inter_stim  # random slope of IV1*IV2 for Stimuli


  # number of levels for IV1 and IV2
  numlevel_IV1 <- length(IV1.levels)
  numlevel_IV2 = length(IV2.levels)
  num_total = numlevel_IV1 * numlevel_IV2 * num_Stim * num_Subj

  # generate factors for all non-dependent variables
  IV1 <- gl(n = numlevel_IV1, k = 1, length = num_total, labels = IV1.levels)
  IV2 <- gl(n = numlevel_IV2, k = numlevel_IV1, length = num_total, labels = IV2.levels)
  Stim <- gl(n = num_Stim, k = numlevel_IV1 * numlevel_IV2, length = num_total, labels = paste0("stim", 1:num_Stim))
  Subj <- gl(n = num_Subj, k = numlevel_IV1 * numlevel_IV2 * num_Stim, length = num_total, labels = paste0("P", 1:num_Subj))

  expdesign <- data.frame(Subject = Subj, Stimuli = Stim, IV2 = IV2, IV1 = IV1)

  # str(expdesign)

  # RT for fixed factors
  modmat.fix <- model.matrix( ~ IV1 * IV2, expdesign)

  # calculate the betas
  intercept = IV1.1_IV2.1
  int_IV1.2 = IV1.2_IV2.1 - intercept
  slp_IV2.2 = (IV1.1_IV2.2 - intercept) / 1  # divided by 1 (dummy coding)
  interaction = (IV1.2_IV2.2 - IV1.1_IV2.2) - int_IV1.2
  coeff.fix <- c(intercept, int_IV1.2, slp_IV2.2, interaction)

  # added random variance (i.e. the residual) to the RT of fixed effect
  RT.fix <- rnorm(n = num_total, mean = modmat.fix %*% coeff.fix, sd = var_residual)

  # Random intercepts of RT for every subject or stimulus
  rnd_int_subj <- rnorm(num_Subj, 0, var_rnd_int_subj)  # Subject
  rnd_int_stim <- rnorm(num_Stim, 0, var_rnd_int_stim)  # Stimulus

  # The sum of the random intercepts of Subjects and Stimuli
  RT.rnd.int <- rnd_int_subj[as.numeric(Subj)] +  rnd_int_stim[as.numeric(Stim)]

  # RT for random slopes
  rnd_slp_IV1_subj <- rnorm(num_Subj, 0, var_rnd_slp_IV1_subj)
  rnd_slp_IV2_subj <- rnorm(num_Subj, 0, var_rnd_slp_IV2_subj)
  rnd_slp_inter_subj <- rnorm(num_Subj, 0, var_rnd_slp_inter_subj)

  rnd_slp_IV1_stim <- rnorm(num_Stim, 0, var_rnd_slp_IV1_stim)
  rnd_slp_IV2_stim <- rnorm(num_Stim, 0, var_rnd_slp_IV2_stim)
  rnd_slp_inter_stim <- rnorm(num_Stim, 0, var_rnd_slp_inter_stim)

  IV1_2 <- paste0("IV1", IV1.levels[2])  # column name
  IV2_2 <- paste0("IV2", IV2.levels[2])  # column name

  IV1.dummy <- modmat.fix[, IV1_2]
  IV2.dummy <- modmat.fix[, IV2_2]
  inter.dummy <- modmat.fix[, paste(IV1_2, IV2_2, sep = ":")]

  RT.rnd.slp <- rnd_slp_IV1_subj[as.numeric(Subj)] * IV1.dummy +
    rnd_slp_IV2_subj[as.numeric(Subj)] * IV2.dummy +
    rnd_slp_inter_subj[as.numeric(Subj)] * inter.dummy +
    rnd_slp_IV1_stim[as.numeric(Stim)] * IV1.dummy +
    rnd_slp_IV2_stim[as.numeric(Stim)] * IV2.dummy +
    rnd_slp_inter_stim[as.numeric(Stim)] * inter.dummy

  # Sum up all the components of RT
  simudata <- data.frame(expdesign, RT = RT.fix + RT.rnd.int + RT.rnd.slp)

  return(simudata)
}


# function to compare the results of fitted model to the pre-set parameters
compare.simu <- function(params, fit, isconfint = FALSE) {
  # author: Haiyang Jin (https://haiyangjin.github.io/)
  # Usage: compare.lmm(params, fit)$fixed
  #        compare.lmm(params, fit)$random
  
  # load library
  library(magrittr)  # pipe
  library(dplyr)
  
  ## load the parameters
  # set the mu for every bin (every condition)
  IV1.1_IV2.1 <- params$IV1.1_IV2.1
  IV1.2_IV2.1 <- params$IV1.2_IV2.1
  IV1.1_IV2.2 <- params$IV1.1_IV2.2
  IV1.2_IV2.2 <- params$IV1.2_IV2.2
  
  # set the variances for lmm (std)
  var_residual <- params$var_residual  # residual
  var_rnd_int_subj <- params$var_rnd_int_subj  # random intercept for Subject
  var_rnd_int_stim <- params$var_rnd_int_stim  # random intercept for Stimuli
  var_rnd_slp_IV1_subj <- params$var_rnd_slp_IV1_subj  # random slope of IV1 for Subject
  var_rnd_slp_IV2_subj <- params$var_rnd_slp_IV2_subj  # random slope of IV2 for Subject
  var_rnd_slp_inter_subj <- params$var_rnd_slp_inter_subj  # random slope of IV1*IV2 for Subject
  var_rnd_slp_IV1_stim <- params$var_rnd_slp_IV1_stim  # random slope of IV1 for Stimuli
  var_rnd_slp_IV2_stim <- params$var_rnd_slp_IV2_stim  # random slope of IV2 for Stimuli
  var_rnd_slp_inter_stim <- params$var_rnd_slp_inter_stim  # random slope of IV1*IV2 for Stimuli
  
  ## Preparation for fixed effect
  # calculate the betas
  intercept = IV1.1_IV2.1
  int_IV1.2 = IV1.2_IV2.1 - intercept
  slp_IV2.2 = (IV1.1_IV2.2 - intercept) / 1  # divided by 1 (dummy coding)
  interaction = (IV1.2_IV2.2 - IV1.1_IV2.2) - int_IV1.2
  coeff.fix <- c(intercept, int_IV1.2, slp_IV2.2, interaction)

  ## Preparation for random effects
  IV1.levels <- params$IV1.levels
  IV2.levels <- params$IV2.levels
  IV1_2 <- paste0("IV1", IV1.levels[2])  # column name
  IV2_2 <- paste0("IV2", IV2.levels[2])  # column name
  
  full.var.df <- data.frame(
    grp = c(rep("Stimuli", 4), rep("Subject", 4), "Residual"),
    var1 = c(rep(c("(Intercept)", IV1_2, IV2_2, paste(IV1_2, IV2_2, sep = ":")), 2), NA)
  )
  
  rnd.var <- {
    as.data.frame(VarCorr(fit)) %>% 
      filter(is.na(var2)) %>% 
      select(grp, var1, sdcor) %>% 
      rename(sd = sdcor)
  }
  
  rnd.var.full <- merge(full.var.df, rnd.var, all = TRUE)
  
  assumed.variance <- c(var_residual,
                        var_rnd_int_stim, var_rnd_slp_IV1_stim, var_rnd_slp_inter_stim, var_rnd_slp_IV2_stim, 
                        var_rnd_int_subj, var_rnd_slp_IV1_subj, var_rnd_slp_inter_subj, var_rnd_slp_IV2_subj
  )
  
  compare.rnd <- data.frame(rnd.var.full, assumed.variance)
  
  
  # add the confident intervals if isconfint == TURE
  if (isconfint) {
    # load library
    library(lme4)
    # confidence intervals for fixed effects
    confint.fix <- as.data.frame(confint.merMod(fit, parm = "beta_", oldNames = FALSE))  # This step takes a little bit longer
    # confidence intervals for random effects
    confint.rnd <- as.data.frame(confint.merMod(fit, parm = "theta_", oldNames = FALSE))  # This step takes too long
    
    # compare fixed
    compare.fixed <- cbind(as.data.frame(fixef(fit)), confint.fix, assumed.coeff = coeff.fix)
    # compare random
    compare.rnd <- cbind(rnd.var.full, confint.rnd, assumed.variance)
    
  } else {
    # compare fixed
    compare.fixed <- cbind(as.data.frame(fixef(fit)), assumed.coeff = coeff.fix)
    # compare random
    compare.rnd <- cbind(rnd.var.full, assumed.variance)
    
  }
  
  compare.lmm <- list(
    fixed = compare.fixed,
    random = compare.rnd
  )
  
  return(compare.lmm)  
}



