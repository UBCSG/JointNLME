bootFun <- function(time, num_patient, censor_value, par = round(md0$fixed_est, 3), disp = round(md0$disp_est, 3), 
                    Sigma, B = 50, i){
  bootEst <- data.frame(matrix(nrow = B, ncol = length(par)))
  bootTi <- data.frame(matrix(nrow = num_patient, ncol = B))
  names(bootEst) <- names(par)
  b = 1
  while (b <= B){
    cat("############### Simulation Iteration: ", i, "; Bootstrap: ", b, "###############","\n")
    data_boot <- create_simdata(time, num_patient, censor_value, true_value_par = par, 
                                disp_value_par = disp, newSigma = Sigma, trt = T)
    data_decay <- data_boot$decay
    data_rebound <- data_boot$rebound
    data_CD4 <- data_boot$CD4
    data_trans <- data_boot$trans
    
    Object_decay <- list(
      modeltype = "nlme",
      response = "y_decay",
      reg_equation = "(eta1 + Lambda.eta1 * tau1i) + eta3 * exp(- (eta2 + Lambda.eta2 * tau2i) * time_decay)",
      distribution = 'normal',
      random_effects = c("tau1i", "tau2i"),
      fixed_param = list(names = c("eta1", "eta2", "eta3"),
                         start_values = c(par$eta1, par$eta2, par$eta3),
                         lower = c(-Inf, -Inf, -Inf),
                         upper = c( Inf, Inf, Inf)),
      disp_param = list(names = c("Lambda.eta1", "Lambda.eta2"),
                        start_values = c(disp$Lambda.eta1, disp$Lambda.eta2),
                        lower = c(0, 0),
                        upper = c(Inf, Inf)),
      sigma = disp$y_decay_sigma,
      left_censoring = list(indicator = "censor_decay",
                            limit_value = censor_value,
                            method = "tobit"))
    
    Object_rebound <- list(
      modeltype = "nlme",
      response = "y_rebound",
      reg_equation = "(beta1 + Lambda.beta1 * b1i) * time_rebound / (time_rebound + exp(beta2 - (beta3 + Lambda.beta3 * b3i) * time_rebound)) + (beta4)",
      distribution = 'normal',
      random_effects = c("b1i", "b3i"),
      fixed_param = list(names = c("beta1","beta2", "beta3", "beta4"),
                         start_values = c(par$beta1, par$beta2, par$beta3, par$beta4),
                         lower = c(-Inf, -Inf, -Inf, -Inf),
                         upper = c( Inf, Inf, Inf, Inf)),
      disp_param = list(names = c("Lambda.beta1",  "Lambda.beta3"),
                        start_values = c(disp$Lambda.beta1, disp$Lambda.beta3),
                        lower = c(0, 0),
                        upper = c(Inf, Inf)),
      sigma = disp$y_rebound_sigma,
      left_censoring = list(indicator = "censor_rebound",
                            limit_value = censor_value,
                            method = "tobit"))
    
    Object_CD4 <- list(
      modeltype = "nlme",
      response = "CD4",
      reg_equation = "(alpha0 + Lambda.alpha0 * a0i) + (alpha1) * time_CD4 + alpha2 * time_CD4 ^ 2",
      distribution = 'normal',
      random_effects = c('a0i'),
      fixed_param = list(names = paste0("alpha", 0:2),
                         start_values = c(par$alpha0, par$alpha1, par$alpha2),
                         lower = c(-Inf, -Inf, -Inf),
                         upper = c( Inf, Inf, Inf)),
      disp_param = list(names = c("Lambda.alpha0"),
                        start_values = c(disp$Lambda.alpha0),
                        lower = c(0),
                        upper = c(Inf)),
      sigma = disp$CD4_sigma)
    
    Object_trans <- list(
      modeltype = "survival",
      response = "time",
      event = "trans",
      reg_equation = "omega1 * tau2i + omega2 * maxCD4",
      distribution = 'log-normal',
      fixed_param = list(names = c("omega1", "omega2"),
                         start_values = c(par$omega1, par$omega2),
                         lower = c(-Inf, -Inf),
                         upper = c(Inf, Inf)),
      sigma = disp$time_sigma,
      truncated = list(lower = "treatment_stop",
                       upper = "first_reboundobs"))
    
    
    memList <- list(Object_decay, Object_rebound, Object_CD4, Object_trans)
    dataList <- list(data_decay, data_rebound, data_CD4, data_trans)
    
    md0 <- try(fit_multi_mems(memList = memList,
                              dataList = dataList,
                              subject_id = "ID",
                              randeff_info = list(distribution = "normal", degree = 0),
                              # cov_method = c("cholesky"),
                              cov_method = "spherical",
                              loglike_tol = 1e-2, par_tol = 1e-2,
                              Silent = T, iterMax = 30, REML = FALSE, naive = FALSE, adjust = TRUE), silent = T) 
    error_mess <- attr(md0, "class")
    if (!("try-error" %in% error_mess)){
      bootEst[b, ] <- round(md0$fixed_est, 3)
      bootTi[, b] <- md0$Ti$Ti0_sim
      b <- b + 1
    }
  }
  res = list(bootEst = bootEst, bootTi = bootTi)
}

