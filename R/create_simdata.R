create_simdata <- function(time, num_patient, censor_value, true_value_par, 
                           disp_value_par, newSigma, trt){
  # Create dataset for viral decay 
  ID <- rep(1:num_patient, each = length(time))
  t <- rep(time, num_patient)
  N <- num_patient * length(time)
  data0 <- as.data.frame(cbind(ID, t))
  
  if (trt == TRUE){
    ntreatment <- 10
    treatment_stop <- time[ntreatment]
    treatment <- rep(c(rep("1", ntreatment), rep("0", length(time) - ntreatment)), num_patient)
    data0 <- as.data.frame(cbind(data0, treatment, treatment_stop))
    data0$treatment_stop <- as.numeric(data0$treatment_stop) 
  } else{
    data0$treatment_stop <- 0
  }
  
  ranef_sim <- rmvnorm(num_patient, mean = rep(0, nrow(newSigma)), sigma = newSigma)
  colnames(ranef_sim) <- c("tau1i_sim", "tau2i_sim", "tau3i_sim", "b1i_sim", "b2i_sim", "b3i_sim", "b4i_sim", "a0i_sim", "a1i_sim", "a2i_sim")
  randef_sim <- as.data.frame(cbind(ID = c(1:num_patient), ranef_sim))
  data <- merge(data0, randef_sim, by = "ID") 
  data <- data %>% 
    mutate(e1 = rnorm(N, 0, disp_value_par$y_decay_sigma)) %>% 
    mutate(e2 = rnorm(N, 0, disp_value_par$y_rebound_sigma)) %>% 
    mutate(e3 = rnorm(N, 0, disp_value_par$CD4_sigma)) %>% 
    mutate(CD4_true = (true_value_par$alpha0 + disp_value_par$Lambda.alpha0 * a0i_sim) + (true_value_par$alpha1 + disp_value_par$Lambda.alpha1 * a1i_sim) * t + (true_value_par$alpha2 + disp_value_par$Lambda.alpha2 * a2i_sim) * t^2) %>% 
    mutate(CD4_obs = CD4_true + e3)
  max_CD4 <- data %>%
    group_by(ID) %>%
    filter(CD4_true == max(CD4_true)) %>%
    dplyr::select(ID, CD4_true) 
  colnames(max_CD4) <- c("ID", "maxCD4")
  
  max_CD4_obs <- data %>% 
    group_by(ID) %>%
    filter(CD4_obs == max(CD4_obs)) %>%
    dplyr::select(ID, CD4_obs) 
  colnames(max_CD4_obs) <- c("ID", "maxCD4_obs")
  
  randef_sim <- merge(randef_sim, max_CD4)
  randef_sim <- merge(randef_sim, max_CD4_obs)
  
  Ti <- randef_sim %>% 
    mutate(e4 = rnorm(num_patient, 0,  disp_value_par$time_sigma)) %>% 
    mutate(logTi0 = true_value_par$omega1 * tau2i_sim + true_value_par$omega2 * maxCD4 + e4) %>% 
    mutate(Ti0 = exp(logTi0)) %>% 
    dplyr::select(c(ID, Ti0, maxCD4, maxCD4_obs))
  
  par(mfrow = c(1, 2))
  boxplot(Ti$Ti0)
  hist(Ti$Ti0, main = "", xlab = "Ti0")
  
  data <- merge(data, Ti, by = "ID") 
  if (trt == TRUE){
    data <- data %>%
      mutate(Ti_true = Ti0 + treatment_stop) 
  } else{
    data <- data %>%
      mutate(Ti_true = Ti0) 
  }
 data <- data %>% 
    mutate(phase_true = ifelse(t <= Ti_true, "decay", "rebound")) %>%
    mutate(viral = ifelse(phase_true == "decay",
                          (true_value_par$eta1 + disp_value_par$Lambda.eta1 * tau1i_sim) + (true_value_par$eta3  + disp_value_par$Lambda.eta3 * tau3i_sim) * exp(- (true_value_par$eta2 + disp_value_par$Lambda.eta2 * tau2i_sim) * t) + e1,
                          (true_value_par$beta1 + disp_value_par$Lambda.beta1 * b1i_sim) * (t - Ti_true)/ ((t - Ti_true) + exp(true_value_par$beta2 + disp_value_par$Lambda.beta2 * b2i_sim - (true_value_par$beta3 + disp_value_par$Lambda.beta3 * b3i_sim) * (t - Ti_true))) + (true_value_par$beta4 + disp_value_par$Lambda.beta4 * b4i_sim) + e2)) %>%
                          # (true_value_par$beta1) * (t - Ti_true)/ ((t - Ti_true) + exp(true_value_par$beta2 - (true_value_par$beta3) * (t - Ti_true))) + (true_value_par$beta4 + disp_value_par$zeta3 * b4i_sim) + e2)) %>%
    mutate(censor = ifelse(viral < censor_value, 1, 0)) %>%
    mutate(viral = ifelse(viral < censor_value, censor_value, viral)) 
  
 last_decayobs <- data %>%
   filter(phase_true == "decay", !censor) %>%   # keep uncensored and above detection
   group_by(ID) %>%
   slice_max(t) %>%
   ungroup() %>% 
   dplyr::select(ID, t)
 colnames(last_decayobs)[2] <- "last_decayobs"
 data <- merge(data, last_decayobs)
 
  first_reboundobs <- data %>%
    filter(phase_true == "rebound", !censor) %>%   # keep uncensored and above detection
    group_by(ID) %>%
    slice_min(t) %>%
    ungroup() %>% 
    dplyr::select(ID, t)
  colnames(first_reboundobs)[2] <- "first_reboundobs"
  data <- merge(data, first_reboundobs)
  
  data_decay <- data %>% 
    dplyr::select(ID, t, viral, phase_true, censor, Ti_true, treatment_stop, last_decayobs, first_reboundobs)
  colnames(data_decay) <- c("ID", "time_decay", "y_decay", "phase_true", "censor_decay", "Ti_true", "treatment_stop", "last_decayobs", "first_reboundobs")
  data_rebound <- data %>% 
    dplyr::select(ID, t, viral, censor, Ti_true, treatment_stop, last_decayobs, first_reboundobs)
  colnames(data_rebound) <- c("ID", "time_rebound", "y_rebound", "censor_rebound", "Ti_true", "treatment_stop", "last_decayobs", "first_reboundobs")
  data_CD4 <- data %>% 
    dplyr::select(ID, t, CD4_obs)
  colnames(data_CD4) <- c("ID", "time_CD4", "CD4")
  
  data_trans <- data %>%
    dplyr::select(ID, Ti0, Ti_true, treatment_stop, last_decayobs, first_reboundobs, maxCD4_obs) %>% 
    distinct(ID, .keep_all = TRUE) 

  if (trt == TRUE){
    data_decay$last_decayobs <- data_rebound$last_decayobs <- data_trans$last_decayobs <- 0
  } else{
    data_decay$treatment_stop <- data_rebound$treatment_stop <- data_trans$treatment_stop <- 0
  }
  
  list(decay = data_decay, rebound = data_rebound, CD4 = data_CD4, trans = data_trans)
}
