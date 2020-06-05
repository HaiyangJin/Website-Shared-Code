# function to simulate the data set for fitting linear mixed model
simudata_lmm <- function (params, usetibble = 0) {
  # This function simulates a dataset for fitting linear mixed mode.
  # A more detailed explanation could be found [here](https://haiyangjin.github.io/portfolio/simulate-data-lmm/)
  # author: Haiyang Jin (https://haiyangjin.github.io/)
  
  # Please make sure library(tibble) is installed if usetibble = 1.
  
  # # HOW to use this function? (Usage)
  # source("path/to/simulate_data_lmm.R")
  # params <- list(
  #   # set the parameters for experiment design
  #   IV1_levels = c("upright", "inverted"),  # (has to be two levels)
  #   IV2_levels = c("intact", "scrambled"),  # (has to be two levels)
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
  # simudata <- simudata_lmm(params = params)
  #
  # library(lme4)
  # lmm <- lmer(Resp ~ IV1 * IV2 + (1 + IV1 * IV2|Subject) + (1 + IV1 * IV2 | Stimuli), data = simudata) # full model
  # summary(lmm)
  
  # load library
  if (usetibble) {library(tibble)}
  
  ## load the parameters from params
  # set the parameters for experiment design
  IV1_levels <- params$IV1_levels
  IV2_levels <- params$IV2_levels
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
  numlevel_IV1 <- length(IV1_levels)
  numlevel_IV2 = length(IV2_levels)
  num_total = numlevel_IV1 * numlevel_IV2 * num_Stim * num_Subj
  
  # generate factors for all non-dependent variables
  IV1 <- gl(n = numlevel_IV1, k = 1, length = num_total, labels = IV1_levels)
  IV2 <- gl(n = numlevel_IV2, k = numlevel_IV1, length = num_total, labels = IV2_levels)
  Stim <- gl(n = num_Stim, k = numlevel_IV1 * numlevel_IV2, length = num_total, labels = paste0("stim", 1:num_Stim))
  Subj <- gl(n = num_Subj, k = numlevel_IV1 * numlevel_IV2 * num_Stim, length = num_total, labels = paste0("P", 1:num_Subj))
  
  if (usetibble) {
    expdesign <- tibble(Subject = Subj, Stimuli = Stim, IV2 = IV2, IV1 = IV1)
  } else {
    expdesign <- data.frame(Subject = Subj, Stimuli = Stim, IV2 = IV2, IV1 = IV1)
  }
  
  # str(expdesign)
  
  # Resp for fixed factors
  modmat_fix <- model.matrix( ~ IV1 * IV2, expdesign)
  
  # calculate the betas
  intercept = IV1.1_IV2.1
  int_IV1.2 = IV1.2_IV2.1 - intercept
  slp_IV2.2 = (IV1.1_IV2.2 - intercept) / 1  # divided by 1 (dummy coding)
  interaction = (IV1.2_IV2.2 - IV1.1_IV2.2) - int_IV1.2
  coeff_fix <- c(intercept, int_IV1.2, slp_IV2.2, interaction)
  
  # added random variance (i.e. the residual) to the Resp of fixed effect
  Resp_fix <- rnorm(n = num_total, mean = modmat_fix %*% coeff_fix, sd = var_residual)
  
  # Random intercepts of Resp for every subject or stimulus
  rnd_int_subj <- rnorm(num_Subj, 0, var_rnd_int_subj)  # Subject
  rnd_int_stim <- rnorm(num_Stim, 0, var_rnd_int_stim)  # Stimulus
  
  # The sum of the random intercepts of Subjects and Stimuli
  RT_rnd_int <- rnd_int_subj[as.numeric(Subj)] +  rnd_int_stim[as.numeric(Stim)]
  
  # RT for random slopes
  rnd_slp_IV1_subj <- rnorm(num_Subj, 0, var_rnd_slp_IV1_subj)
  rnd_slp_IV2_subj <- rnorm(num_Subj, 0, var_rnd_slp_IV2_subj)
  rnd_slp_inter_subj <- rnorm(num_Subj, 0, var_rnd_slp_inter_subj)
  
  rnd_slp_IV1_stim <- rnorm(num_Stim, 0, var_rnd_slp_IV1_stim)
  rnd_slp_IV2_stim <- rnorm(num_Stim, 0, var_rnd_slp_IV2_stim)
  rnd_slp_inter_stim <- rnorm(num_Stim, 0, var_rnd_slp_inter_stim)
  
  IV1_2 <- paste0("IV1", IV1_levels[2])  # column name
  IV2_2 <- paste0("IV2", IV2_levels[2])  # column name
  
  IV1_dummy <- modmat_fix[, IV1_2]
  IV2_dummy <- modmat_fix[, IV2_2]
  inter_dummy <- modmat_fix[, paste(IV1_2, IV2_2, sep = ":")]
  
  RT_rnd_slp <- rnd_slp_IV1_subj[as.numeric(Subj)] * IV1_dummy +
    rnd_slp_IV2_subj[as.numeric(Subj)] * IV2_dummy +
    rnd_slp_inter_subj[as.numeric(Subj)] * inter_dummy +
    rnd_slp_IV1_stim[as.numeric(Stim)] * IV1_dummy +
    rnd_slp_IV2_stim[as.numeric(Stim)] * IV2_dummy +
    rnd_slp_inter_stim[as.numeric(Stim)] * inter_dummy
  
  # Sum up all the components of RT
  if (usetibble) {
    simudata <- add_column(expdesign, RT = RT_fix + RT_rnd_int + RT_rnd_slp)
  } else {
    simudata <- data.frame(expdesign, RT = RT_fix + RT_rnd_int + RT_rnd_slp)
  }
  return(simudata)
}


# function to compare the results of fitted model to the pre-set parameters
compare_simu <- function(params, fit, isconfint = FALSE) {
  # author: Haiyang Jin (https://haiyangjin.github.io/)
  # Usage: compare_simu(params, fit)$fixed
  #        compare_simu(params, fit)$random
  
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
  coeff_fix <- c(intercept, int_IV1.2, slp_IV2.2, interaction)
  
  ## Preparation for random effects
  IV1_levels <- params$IV1_levels
  IV2_levels <- params$IV2_levels
  IV1_2 <- paste0("IV1", IV1_levels[2])  # column name
  IV2_2 <- paste0("IV2", IV2_levels[2])  # column name
  
  full_var_df <- data.frame(
    grp = c(rep("Stimuli", 4), rep("Subject", 4), "Residual"),
    var1 = c(rep(c("(Intercept)", IV1_2, IV2_2, paste(IV1_2, IV2_2, sep = ":")), 2), NA)
  )
  
  rnd_var <- {
    as.data.frame(VarCorr(fit)) %>% 
      filter(is.na(var2)) %>% 
      select(grp, var1, sdcor) %>% 
      rename(sd = sdcor)
  }
  
  rnd_var_full <- merge(full_var_df, rnd_var, all = TRUE)
  
  assumed_variance <- c(var_residual,
                        var_rnd_int_stim, var_rnd_slp_IV1_stim, var_rnd_slp_inter_stim, var_rnd_slp_IV2_stim, 
                        var_rnd_int_subj, var_rnd_slp_IV1_subj, var_rnd_slp_inter_subj, var_rnd_slp_IV2_subj
  )
  
  compare_rnd <- data.frame(rnd_var_full, assumed_variance)
  
  
  # add the confident intervals if isconfint == TURE
  if (isconfint) {
    # load library
    library(lme4)
    # confidence intervals for fixed effects
    confint_fix <- as.data.frame(confint.merMod(fit, parm = "beta_", oldNames = FALSE))  # This step takes a little bit longer
    # confidence intervals for random effects
    confint_rnd <- as.data.frame(confint.merMod(fit, parm = "theta_", oldNames = FALSE))  # This step takes too long
    
    # compare fixed
    compare_fixed <- cbind(as.data.frame(fixef(fit)), confint_fix, assumed_coeff = coeff_fix)
    # compare random
    compare_rnd <- cbind(rnd_var_full, confint_rnd, assumed_variance)
    
  } else {
    # compare fixed
    compare_fixed <- cbind(as.data.frame(fixef(fit)), assumed_coeff = coeff_fix)
    # compare random
    compare_rnd <- cbind(rnd_var_full, assumed_variance)
    
  }
  
  compare_lmm <- list(
    fixed = compare_fixed,
    random = compare_rnd
  )
  
  return(compare_lmm)  
}



