create_simdata_new <- function(time, num_patient, censor_value, true_value_par, 
                           disp_value_all, Sigma){
  # Create dataset for viral decay 
  ID <- rep(1:num_patient, each = length(time))
  t <- rep(time, num_patient)
  N <- num_patient * length(time)
  data0 <- as.data.frame(cbind(ID, t))
  
  ranef_sim <- rmvnorm(num_patient, mean = rep(0, nrow(Sigma)), sigma = Sigma)
  colnames(ranef_sim) <- c("tau1i_sim", "tau2i_sim", "tau3i_sim", "b1i_sim", "b2i_sim", "b3i_sim", "b4i_sim", "a0i_sim", "a1i_sim", "a2i_sim")
  randef_sim <- as.data.frame(cbind(ID = c(1:num_patient), ranef_sim))
  data <- merge(data0, randef_sim, by = "ID") 
  data <- data %>% 
    mutate(e1 = rnorm(N, 0, disp_value_all$y_decay_sigma)) %>% 
    mutate(e2 = rnorm(N, 0, disp_value_all$y_rebound_sigma)) %>% 
    mutate(e3 = rnorm(N, 0, disp_value_all$CD4_sigma)) %>% 
    mutate(CD4_true = (true_value_par$alpha0 + disp_value_all$Lambda.alpha0 * a0i_sim) + (true_value_par$alpha1 + disp_value_all$Lambda.alpha1 * a1i_sim) * t + (true_value_par$alpha2 + disp_value_all$Lambda.alpha2 * a2i_sim) * t^2) %>% 
    mutate(CD4_obs = CD4_true + e3)
  max_CD4 <- data %>%
    group_by(ID) %>%
    filter(CD4_true == max(CD4_true)) %>%
    dplyr::select(ID, CD4_true) 
    
  colnames(max_CD4) <- c("ID", "maxCD4")
  randef_sim <- merge(randef_sim, max_CD4)
  
  Ti <- randef_sim %>% 
    mutate(e4 = rnorm(num_patient, 0,  disp_value_all$time_sigma)) %>% 
    mutate(logTi_true = true_value_par$omega1 * tau2i_sim + true_value_par$omega2 * maxCD4 + e4) %>% 
    mutate(Ti_true = exp(logTi_true)) %>% 
    dplyr::select(c(ID, Ti_true, maxCD4)) 
  
  # Ti <- randef_sim %>% 
  #   group_by(ID) %>% 
  #   mutate(Ti_true = rlnorm(1,  true_value_par$omega1 * tau2i_sim + true_value_par$omega2 * maxCD4, disp_value_all$time_sigma)) %>% 
  #   dplyr::select(c(ID, Ti_true, maxCD4))
  
  par(mfrow = c(1, 2))
  boxplot(Ti$Ti_true)
  hist(Ti$Ti_true, main = "", xlab = "Ti")
  
  data <- merge(data, Ti, by = "ID") 
  data <- data %>%
    mutate(phase = ifelse(t <= Ti_true, "decay", "rebound")) %>%
    mutate(viral = ifelse(phase == "decay",
                          (true_value_par$eta1 + disp_value_all$Lambda.eta1 * tau1i_sim) + (true_value_par$eta3  + disp_value_all$Lambda.eta3 * tau3i_sim) * exp(- (true_value_par$eta2 + disp_value_all$Lambda.eta2 * tau2i_sim) * t) + e1,
                          (true_value_par$beta1 + disp_value_all$Lambda.beta1 * b1i_sim) * (t - Ti_true)/ ((t - Ti_true) + exp(true_value_par$beta2 + disp_value_all$Lambda.beta2 * b2i_sim - (true_value_par$beta3 + disp_value_all$Lambda.beta3 * b3i_sim) * (t - Ti_true))) + (true_value_par$beta4 + disp_value_all$Lambda.beta4 * b4i_sim) + e2)) %>%
                          # (true_value_par$beta1) * (t - Ti_true)/ ((t - Ti_true) + exp(true_value_par$beta2 - (true_value_par$beta3) * (t - Ti_true))) + (true_value_par$beta4 + disp_value_all$zeta3 * b4i_sim) + e2)) %>%
    mutate(censor = ifelse(viral < censor_value, 1, 0)) %>%
    mutate(viral = ifelse(viral < censor_value, censor_value, viral))
  
  last_decayobs <- data %>%
    filter(phase == "decay", !censor) %>%   # keep uncensored and above detection
    group_by(ID) %>%
    slice_max(t) %>%
    ungroup() %>% 
    dplyr::select(ID, t)
  colnames(last_decayobs)[2] <- "last_decayobs"
  data <- merge(data, last_decayobs)
  
  first_reboundobs <- data %>%
    filter(phase == "rebound", !censor) %>%   # keep uncensored and above detection
    group_by(ID) %>%
    slice_min(t) %>%
    ungroup() %>% 
    dplyr::select(ID, t)
  colnames(first_reboundobs)[2] <- "first_reboundobs"
  data <- merge(data, first_reboundobs)
  
  data_decay <- data %>% 
    dplyr::select(ID, t, phase, viral, censor, Ti_true, last_decayobs, first_reboundobs)
  colnames(data_decay) <- c("ID", "time_decay", "phase", "y_decay", "censor_decay", "Ti_true", "last_decayobs", "first_reboundobs")
  data_rebound <- data %>% 
    dplyr::select(ID, t, phase, viral, censor, Ti_true, last_decayobs, first_reboundobs)
  colnames(data_rebound) <- c("ID", "time_rebound", "phase", "y_rebound", "censor_rebound", "Ti_true", "last_decayobs", "first_reboundobs")
  data_CD4 <- data %>% 
    dplyr::select(ID, t, CD4_obs)
  colnames(data_CD4) <- c("ID", "time_CD4", "CD4")
  
  data_trans <- data %>%
    dplyr::select(ID, Ti_true, last_decayobs, first_reboundobs) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    mutate(trans = 1)

  list(decay = data_decay, rebound = data_rebound, CD4 = data_CD4, trans = data_trans)
}
