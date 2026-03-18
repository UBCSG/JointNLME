library(ggplot2)
library(magrittr)
library(lme4)
library(readxl)
library(dplyr)
library(tidyverse)
library(nlme)
library(Deriv)
library(matrixcalc)
library(optimx)
library(MASS)
library(survival)
library(ReIns)
library(mvtnorm)
library(icenReg)
library(MBESS)
library(xtable)
library(EnvStats)
library(gdata)

setwd("R")
file.sources = list.files(pattern = "*.R")
sapply(file.sources, source, .GlobalEnv)

setwd(here::here())

# Some true values
num_sim <- 100  # number of simulation
num_patient <- 100  # number of patients
censor_value <- log10(40)  # censor value 

# Time points 
# time <- c(0.03, 0.27,0.47, 0.93, 1.25, 1.43, 1.87, 2.17, 2.25, 2.56, 2.80, 3.00, 3.15, 3.34, 3.40, 3.68,
#           3.87, 4.00, 4.2, 4.37, 4.65, 4.78, 4.94, 5.04, 5.20, 5.45, 6.00, 6.10, 6.31, 6.53, 6.87, 7.10,
#           7.28, 7.45, 7.76, 8.1, 8.42, 8.8, 9.1)
time <- c(0.03, 0.27,0.47, 0.93, 1.25, 1.43, 1.87, 2.25, 2.80, 3.15, 3.40, 
          3.87,  4.2, 4.65, 4.94, 5.20, 5.45, 6.00, 6.31, 6.53, 6.87, 
          7.28, 7.76, 8.1, 8.42, 9.1)

# Viral decay model parameter
eta1 = 0.8
eta2 = 0.7
eta3 = 3
sigma1 <- 0.06
Lambda.eta <- c(0.1, 0.05, 0)
# Viral rebound model parameter 
beta1 = 1.8
beta2 = 1
beta3 = 3
beta4 = 1.2
sigma2 <- 0.1
Lambda.beta <- c(0.1, 0, 0.1, 0)
# CD4 model parameter
alpha0 = 3
alpha1 = 0.2
alpha2 = -0.03
sigma3 <- 0.4
Lambda.alpha <- c(0.1, 0, 0)
# Transition point model parameter
omega1 <- -0.2
omega2 <- 0.1
sigma4 <- 0.1

names(Lambda.eta) <- paste0("Lambda.eta", 1:length(Lambda.eta))
names(Lambda.beta) <- paste0("Lambda.beta", 1:length(Lambda.beta))
names(Lambda.alpha) <- paste0("Lambda.alpha", 0:(length(Lambda.alpha)-1))
# Combine all into one named vector
Lambda.all <- c(Lambda.eta, Lambda.beta, Lambda.alpha)
# Filter out values that are not zero
Lambda.nonzero <- Lambda.all[Lambda.all != 0]
Lambda.zero <- Lambda.all[Lambda.all == 0]
pos.Lambda.nonzero <- which(Lambda.all != 0)

n_randef <- length(Lambda.nonzero)
Sigma = as.matrix(diag(rep(1, n_randef))) 
Sigma[2, 3] <- Sigma[3, 2] <- 0.8
newSigma <- matrix(0, nrow = 10, ncol = 10)
newSigma[pos.Lambda.nonzero, pos.Lambda.nonzero] <- Sigma
diag(newSigma) <- 1

# Results Table
upper_Sigma <- Sigma
upper_Sigma[lower.tri(upper_Sigma, diag = TRUE)] <- 0
nonzero_indices <- which(upper_Sigma != 0, arr.ind = TRUE)
Sigma_names <- apply(nonzero_indices, 1, function(idx) paste0("Sigma", idx[1], idx[2]))

par_names <- c("eta1", "eta2", "eta3", "beta1", "beta2", "beta3", "beta4", "alpha0", "alpha1", "alpha2", "omega1", "omega2")
par_names_naive <- c("eta1", "eta2", "eta3", "beta1", "beta2", "beta3", "beta4", "omega1", "omega2")
disppar_names <- c(Sigma_names, names(Lambda.nonzero), "y_decay_sigma", "y_rebound_sigma", "CD4_sigma", "time_sigma")
disppar_names_naive <- disppar_names[-c(6, 9)]

true_value <- c(eta1, eta2, eta3, beta1, beta2, beta3, beta4, alpha0, alpha1, alpha2, omega1, omega2)
par_length <- length(true_value)
disp_value <- c(Sigma[nonzero_indices], Lambda.nonzero, sigma1, sigma2, sigma3, sigma4)
disppar_length <- length(disp_value)
true_value_par <- as.list(true_value)
names(true_value_par) <- par_names
disp_value_par <- as.list(c(disp_value, Lambda.zero))
names(disp_value_par)[1: disppar_length] <- disppar_names
disp_value_naive <- disp_value[-c(6,9)]

table_par <- table_parsd <- boot_sd <- data.frame(matrix(nrow = num_sim, ncol = par_length))
table_disp <- data.frame(matrix(nrow = num_sim, ncol = disppar_length))
colnames(table_par) <- colnames(table_parsd) <- colnames(boot_sd) <- c(par_names)
colnames(table_disp) <- disppar_names

true_value_naive <- c(eta1, eta2, eta3, beta1, beta2, beta3, beta4, omega1, omega2)
par_length_naive <- length(true_value_naive)
disppar_length_naive <- length(disppar_names_naive)
table_par_naive <- table_parsd_naive <- boot_sd_naive <- data.frame(matrix(nrow = num_sim, ncol = par_length_naive))
colnames(table_par_naive) <- colnames(table_parsd_naive) <- colnames(boot_sd_naive) <- c(par_names_naive)
table_disp_naive <- data.frame(matrix(nrow = num_sim, ncol = disppar_length_naive))
colnames(table_disp_naive) <- disppar_names_naive

table_TiCI <- data.frame(matrix(nrow = num_sim, ncol = 3))
colnames(table_TiCI) <- c("New", "Naive", "Poly")
table_Ti_list <- list()
  
table_par <- readRDS("SimResults/table_par.RDS")
table_disp <- readRDS("SimResults/table_disp.RDS")
table_parsd <- readRDS("SimResults/table_parsd.RDS")
table_par_naive <- readRDS("SimResults/table_par_naive.RDS")
table_disp_naive <- readRDS("SimResults/table_disp_naive.RDS")
table_parsd_naive <- readRDS("SimResults/table_parsd_naive.RDS")
table_Ti_list <- readRDS("SimResults/table_Ti_list.RDS")
table_TiCI <- readRDS("SimResults/table_TiCI.RDS")


i = 1
NC <- 0

while (i <= num_sim){
  print(glue::glue("i = ",i, " New method"))
  
  data_sim <- create_simdata(time, num_patient, censor_value, true_value_par, 
                             disp_value_par, newSigma, trt = T)
  
  
  data_decay <- data_sim$decay
  data_rebound <- data_sim$rebound
  data_CD4 <- data_sim$CD4
  data_trans <- data_sim$trans
  
  mean(data_decay$censor)
  # ggplot(data_decay[data_decay$ID %in% c(1:10),], aes(x=time_decay, y = y_decay)) +
  #   geom_point(aes(fill=factor(censor_decay)),size=2, shape=21, stroke=0) +
  #   geom_line(aes(group=ID, color = factor(phase_true))) +
  #   scale_x_continuous("Time") +
  #   scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  #   scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
  #   scale_color_manual(values=c("black", "steelblue"), labels = c("After treatment interruption", "During treatment"))+
  #   labs(color = "Treatment", fill = "Data type")+ggtitle("Viral load trajectories of 10 randomly selected subjects") +
  #   theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
  # 
  subdata <- data_decay[data_decay$ID %in% c(1),]
  tq1_trt <- unique(subdata$treatment_stop)
  tq1_lastdecayobs <- unique(subdata$last_decayobs)
  tq1 <- max(tq1_lastdecayobs, tq1_trt)
  tq2 <- unique(subdata$first_reboundobs)
  Ti <- unique(data_trans[data_trans$ID %in% c(1),]$Ti_true)
  # ggplot(subdata, aes(x=time_decay, y = y_decay)) +
  #   geom_point(aes(fill=factor(censor_decay)),size=2, shape=21, stroke=0) +
  #   geom_line() +
  #   scale_x_continuous("Time") +
  #   geom_vline(xintercept = tq2, linetype = 2) +
  #   annotate("text", x = 1.13 * tq2, y = 3.5, label= "t[iq2]",
  #            parse=TRUE, hjust = 1.1) +
  #   geom_vline(xintercept = tq1, linetype = 2) +
  #   annotate("text", x = 0.92 * tq1, y = 3.5, label= "t[iq1]",
  #            parse=TRUE, hjust = 1.1) +
  #   geom_vline(xintercept = Ti, linetype = 2) +
  #   annotate("text", x = 1.1 * Ti, y = 3.5, label= "T[i]",
  #            parse=TRUE, hjust = 1.1) +
  #   scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  #   scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
  #   labs(color = "Phase", fill = "Data type")+ggtitle("Viral load trajectory of a typical individual") +
  #   theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
  # 
  # ggplot(data_CD4[data_CD4$ID %in%c(1:10),], aes(x = time_CD4, y = CD4)) +
  #   geom_line(aes(group=ID)) +
  #   scale_x_continuous("Time") +
  #   scale_y_continuous(bquote("CD4 values"))+
  #   scale_color_manual(values=c("steelblue","black"),labels = c("Following ART interruption", "During ART"))+
  #   labs(color = "ART status")+ggtitle("CD4 plot for 10 randomly selected individuals")+
  #   theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))

  
  #### =============== New method ==================
  Object_decay <- list(
    modeltype = "nlme",
    response = "y_decay",
    reg_equation = "(eta1 + Lambda.eta1 * tau1i) + eta3 * exp(- (eta2 + Lambda.eta2 * tau2i) * time_decay)",
    distribution = 'normal',
    random_effects = c("tau1i", "tau2i"),
    fixed_param = list(names = c("eta1", "eta2", "eta3"),
                       start_values = c(true_value_par$eta1, true_value_par$eta2, true_value_par$eta3),
                       lower = c(-Inf, -Inf, -Inf),
                       upper = c( Inf, Inf, Inf)),
    disp_param = list(names = c("Lambda.eta1", "Lambda.eta2"),
                      start_values = c(Lambda.eta[1], Lambda.eta[2]),
                      lower = c(0, 0),
                      upper = c(Inf, Inf)),
    sigma = disp_value_par$y_decay_sigma,
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
                       start_values = c(true_value_par$beta1, true_value_par$beta2, true_value_par$beta3, true_value_par$beta4),
                       lower = c(-Inf, -Inf, -Inf, -Inf),
                       upper = c( Inf, Inf, Inf, Inf)),
    disp_param = list(names = c("Lambda.beta1",  "Lambda.beta3"),
                      start_values = c(Lambda.beta[1], Lambda.beta[3]),
                      lower = c(0, 0),
                      upper = c(Inf, Inf)),
    sigma = disp_value_par$y_rebound_sigma,
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
                       start_values = c(true_value_par$alpha0, true_value_par$alpha1, true_value_par$alpha2),
                       lower = c(-Inf, -Inf, -Inf),
                       upper = c( Inf, Inf, Inf)),
    disp_param = list(names = c("Lambda.alpha0"),
                      start_values = c(Lambda.alpha[1]),
                      lower = c(0),
                      upper = c(Inf)),
    sigma = disp_value_par$CD4_sigma)
  
  Object_trans <- list(
    modeltype = "survival",
    response = "time",
    event = "trans",
    reg_equation = "omega1 * tau2i + omega2 * maxCD4",
    distribution = 'log-normal',
    fixed_param = list(names = c("omega1", "omega2"),
                       start_values = c(true_value_par$omega1, true_value_par$omega2),
                       lower = c(-Inf, -Inf),
                       upper = c(Inf, Inf)),
    sigma = disp_value_par$time_sigma,
    truncated = list(lower = "treatment_stop",
                     upper = "first_reboundobs"))
  
  memList <- list(Object_decay, Object_rebound, Object_CD4, Object_trans)
  dataList <- list(data_decay, data_rebound, data_CD4, data_trans)
  
  md0 <- try(fit_multi_mems(memList = memList,
                            dataList = dataList,
                            subject_id = "ID",
                            # randeff_info = list(distribution = "normal", degree = 0, correlation = "independent"),
                            randeff_info = list(distribution = "normal", degree = 0),
                            # cov_method = c("cholesky"),
                            cov_method = "spherical",
                            loglike_tol = 1e-3, par_tol = 1e-3,
                            Silent = T, iterMax = 30, REML = FALSE, naive = FALSE, adjust = TRUE), silent = T) 
  # (true_value > md0$fixed_est - 1.96 * md0$fixed_sd & true_value < md0$fixed_est + 1.96 * md0$fixed_sd)
  error_mess <- attr(md0, "class")
  if ("try-error" %in% error_mess){
    NC <- NC + 1
    next
  } else  if (any(is.na(md0$fixed_sd))){
    NC <- NC + 1
    next
  } else{
    randef <- md0$Bi
    fixed_est <- as.list(md0$fixed_est)
    disp_est <- as.list(c(md0$disp_est, Lambda.zero))
    Ti <- merge(data_CD4, randef) %>%
      mutate(fitted_CD4 = fixed_est$alpha0 + disp_est$Lambda.alpha0 * a0i + fixed_est$alpha1 * time_CD4 + fixed_est$alpha2 * time_CD4 ^ 2) %>%
      group_by(ID) %>%
      mutate(max_CD4 = max(fitted_CD4)) %>%
      dplyr::select(-c(time_CD4, CD4, fitted_CD4)) %>%
      distinct() %>%
      mutate(Ti0_est = exp(fixed_est$omega1 * tau2i + fixed_est$omega2 * max_CD4))
    Ti <- merge(Ti, md0$Ti)
    TiCI <- merge(Ti, data_trans) %>%
      # mutate(lowerCI = exp(log(Ti0_sim) - 1.96 * disp_est$time_sigma)) %>% 
      # mutate(upperCI = exp(log(Ti0_sim) + 1.96 * disp_est$time_sigma)) %>% 
      mutate(lowerCI = Ti0_sim * exp(- 1.96 * disp_est$time_sigma)) %>% 
      mutate(upperCI = Ti0_sim * exp(1.96 * disp_est$time_sigma)) %>% 
      mutate(CI = Ti0 >= lowerCI & Ti0 <= upperCI) %>%
      dplyr::select(Ti0, Ti_true, Ti0_est, Ti0_sim, lowerCI, upperCI, CI)
    
    table_par[i, ] <- round(md0$fixed_est, 4)
    table_disp[i, ] <- c(solve(md0$invSIGMA)[2, 3], round(md0$disp_est, 4))
    table_parsd[i, ] <- round(as.numeric(md0$fixed_sd), 4)
    table_TiCI[i, 1] <- mean(TiCI$CI)
  }
  
  #### =============== Naive ==================
  # print(glue::glue("i = ",i, " Naive method"))
  # Object_trans_naive <- list(
  #   modeltype = "survival",
  #   response = "time",
  #   event = "trans",
  #   reg_equation = "omega1 * tau2i + omega2 * maxCD4_obs",
  #   distribution = 'log-normal',
  #   fixed_param = list(names = c("omega1", "omega2"),
  #                      start_values = c(true_value_par$omega1, true_value_par$omega2),
  #                      lower = c(-Inf, -Inf),
  #                      upper = c(Inf, Inf)),
  #   sigma = disp_value_par$time_sigma,
  #   truncated = list(lower = "treatment_stop",
  #                    upper = "first_reboundobs"))
  # memList_naive <- list(Object_decay, Object_rebound, Object_trans_naive)
  # dataList_naive <- list(data_decay, data_rebound, data_trans)
  # 
  # md0_naive <- try(fit_multi_mems(memList = memList_naive,
  #                           dataList = dataList_naive,
  #                           subject_id = "ID",
  #                           # randeff_info = list(distribution = "normal", degree = 0, correlation = "independent"),
  #                           randeff_info = list(distribution = "normal", degree = 0),
  #                           # cov_method = c("cholesky"),
  #                           cov_method = "spherical",
  #                           loglike_tol = 1e-3, par_tol = 1e-3,
  #                           Silent = T, iterMax = 30, REML = FALSE, naive = TRUE, adjust = TRUE), silent = T) 
  # # (true_value_naive > md0_naive$fixed_est - 1.96 * md0_naive$fixed_sd & true_value_naive < md0_naive$fixed_est + 1.96 * md0_naive$fixed_sd)
  # error_mess <- attr(md0_naive, "class")
  # if ("try-error" %in% error_mess){
  #   next
  # } else{
  #   randef_naive <- md0_naive$Bi
  #   fixed_est_naive <- as.list(md0_naive$fixed_est)
  #   disp_est_naive <- as.list(c(md0_naive$disp_est, Lambda.zero))
  #   # Ti_naive <- merge(data_trans, md0_naive$Bi) %>%
  #   #   mutate(Ti_est_naive = exp(fixed_est_naive$omega1 * tau2i + fixed_est_naive$omega2 * maxCD4_obs)) %>%
  #   #   mutate(lowerCI_naive = qlnormTrunc(0.025, meanlog = log(Ti_est_naive), sdlog = disp_est_naive$time_sigma, min = last_decayobs, max = first_reboundobs - treatment_stop)) %>%
  #   #   mutate(upperCI_naive = qlnormTrunc(0.975, meanlog = log(Ti_est_naive), sdlog = disp_est_naive$time_sigma, min = last_decayobs, max = first_reboundobs - treatment_stop)) %>%
  #   #   mutate(CI_naive = Ti0 >= lowerCI_naive & Ti0 <= upperCI_naive) %>%
  #   #   dplyr::select(Ti_est_naive, lowerCI_naive, upperCI_naive, CI_naive)
  #   Ti_naive <- merge(data_trans, md0_naive$Bi) %>%
  #     mutate(Ti0_est_naive = exp(fixed_est_naive$omega1 * tau2i + fixed_est_naive$omega2 * maxCD4_obs)) %>% 
  #     mutate(lowerCI_naive = exp(log(Ti0_est_naive) - 1.96 * disp_est_naive$time_sigma)) %>% 
  #     mutate(upperCI_naive = exp(log(Ti0_est_naive) + 1.96 * disp_est_naive$time_sigma)) %>% 
  #     mutate(CI_naive = Ti0 >= lowerCI_naive & Ti0 <= upperCI_naive) %>%
  #     dplyr::select(Ti0_est_naive, lowerCI_naive, upperCI_naive, CI_naive)
  #   
  #   table_par_naive[i, ] <- round(md0_naive$fixed_est, 4)
  #   table_disp_naive[i, ] <- c(solve(md0_naive$invSIGMA)[2, 3], round(md0_naive$disp_est, 4))
  #   table_parsd_naive[i, ] <- round(as.numeric(md0_naive$fixed_sd), 4)
  #   table_TiCI[i, 2] <- mean(Ti_naive$CI_naive)
  # }
  
  #### =============== Polynomial ==================
  data_poly <- data_decay %>%
    arrange(ID, time_decay) %>%
    group_by(ID) %>%
    mutate(
      rebound_count = cumsum(phase_true == "rebound" & censor_decay == 0)
    ) %>%
    filter(rebound_count <= 8 | phase_true == "decay") %>% 
    mutate(y_poly = ifelse(censor_decay == 1, log10(20), y_decay))
  data_poly <- groupedData(y_poly ~ time_decay | ID, data = data_poly)
  model_poly <- lme(y_poly ~ time_decay + I(time_decay^2), data = data_poly, 
                    random = ~ time_decay - 1 | ID, method = "ML")
  # Fixed effects
  a_fixed <- fixed.effects(model_poly)["(Intercept)"]
  b_fixed <- fixed.effects(model_poly)["time_decay"]
  c_fixed <- fixed.effects(model_poly)["I(time_decay^2)"]
  # Extract random effects
  ranefs <- ranef(model_poly)  # Data frame with rownames = patient IDs
  bi <- rownames_to_column(ranefs, "ID")
  # Extract variance of random slope for time_decay
  # (Assumes time_decay is the only random slope and has no correlation with other REs)
  var_bi <- as.numeric(VarCorr(model_poly)["time_decay", "StdDev"])^2
  # Fixed effect variance of c
  var_c <- vcov(model_poly)["I(time_decay^2)", "I(time_decay^2)"]
  # Covariance between b and c (from fixed effects only — no correlation assumed between REs and FE)
  cov_bc <- vcov(model_poly)["time_decay", "I(time_decay^2)"]
  # Compute Ti_hat and CI for each patient
  Ti_poly_table <- data.frame(ID = rownames(ranefs),
                              Ti_poly = NA,
                              lowerCI_poly = NA,
                              upperCI_poly = NA) 
  Ti_poly_table$ID <- as.numeric(Ti_poly_table$ID)
  Ti_poly_table <- Ti_poly_table %>% arrange(ID)
  for (id in rownames(ranefs)) {
    # Patient-specific b_i
    b_i <- b_fixed + bi[which(bi$ID == id), ]$time_decay
    # Ti_hat for this patient
    Ti_poly <- -b_i / (2 * c_fixed)
    # Delta method: gradient is [-1/(2c), b/(2c^2)]
    dTi_db <- -1 / (2 * c_fixed)
    dTi_dc <- b_i / (2 * c_fixed^2)
    grad <- c(dTi_db, dTi_dc)
    cov_mat <- vcov(model_poly)[c("time_decay", "I(time_decay^2)"),
                                c("time_decay", "I(time_decay^2)")]
    var_Ti <- t(grad) %*% cov_mat %*% grad
    se_Ti <- sqrt(as.numeric(var_Ti))
    # Variance of Ti_hat (includes fixed c and random b_i uncertainty)
    # var_Ti <- dTi_db^2 * var_bi + dTi_dc^2 * var_c + 2 * dTi_db * dTi_dc * cov_bc
    # se_Ti <- sqrt(var_Ti)
    # 95% CI
    Ti_poly_table[Ti_poly_table$ID == id, "Ti_poly"] <- Ti_poly
    Ti_poly_table[Ti_poly_table$ID == id, "lowerCI_poly"] <- Ti_poly - 1.96 * se_Ti
    Ti_poly_table[Ti_poly_table$ID == id, "upperCI_poly"] <- Ti_poly + 1.96 * se_Ti
  }
  table_Ti_list[[i]] <- cbind(TiCI, Ti_naive, Ti_poly_table[, c(2:4)]) %>% 
    mutate(CI_poly = Ti_true >= lowerCI_poly & Ti_true <= upperCI_poly)
  table_TiCI[i, 3] <- mean(table_Ti_list[[i]]$CI_poly)
  
  saveRDS(table_par, "SimResults IV/table_par.RDS")
  saveRDS(table_disp, "SimResults IV/table_disp.RDS")
  saveRDS(table_parsd, "SimResults IV/table_parsd.RDS")
  saveRDS(table_par_naive, "SimResults IV/table_par_naive.RDS")
  saveRDS(table_disp_naive, "SimResults IV/table_disp_naive.RDS")
  saveRDS(table_parsd_naive, "SimResults IV/table_parsd_naive.RDS")
  saveRDS(table_Ti_list, "SimResults IV/table_Ti_list.RDS")
  saveRDS(table_TiCI, "SimResults IV/table_TiCI.RDS")
  i <- i + 1
}

q = 1
randef <- as.list(md0$Bi %>% filter(ID == q))
fixed_est <- as.list(md0$fixed_est)
subTidata <- data_trans %>% filter(ID == q)
Ti0_est <- TiCI[q, ]$Ti0_sim
Ti_est <- Ti0_est + tq1_trt
Ti_est_lower <- TiCI[q, ]$lowerCI + tq1_trt
Ti_est_upper <- TiCI[q, ]$upperCI + tq1_trt
if (Ti_est_lower < tq1){
  Ti_est_lower <- tq1
}
if (Ti_est_upper > tq2){
  Ti_est_upper <- tq2
}
# Ti_est_poly <- Ti_poly_table[q, ]$Ti_poly
trt_stop <- unique(subdata$treatment_stop)
(Ti_bias = subTidata$Ti_true - Ti_est)
# (Ti_poly_bias = subTidata$Ti_true - Ti_est_poly)
# b_i <- b_fixed + bi[which(bi$ID == q), ]$time_decay
t_end <- max(subdata$time_decay)
subdata <- subdata %>% 
  mutate(y_decay_halfDL = ifelse(censor_decay == 1, log10(1/2*10^(y_decay)), y_decay))
ggplot(subdata, aes(x = time_decay, y = y_decay_halfDL)) +
  geom_point(aes(shape = factor(censor_decay)), size = 2) +
  scale_shape_manual(name = "Censored", values=c(19, 1), labels = c("No", "Yes")) +
  geom_function(fun = function(x) (fixed_est$eta1 + disp_est$Lambda.eta1 * randef$tau1i) + fixed_est$eta3 * exp(- (fixed_est$eta2 + disp_est$Lambda.eta2 * randef$tau2i) * x), 
                xlim = c(0, 2.2)) + 
  geom_function(fun = function(x) ((fixed_est$beta1 + disp_est$Lambda.beta1 * randef$b1i) * (x - Ti_est) / ((x - Ti_est) + exp(fixed_est$beta2 - (fixed_est$beta3 + disp_est$Lambda.beta3 * randef$b3i) * (x - Ti_est))) + fixed_est$beta4),
                xlim = c(Ti_est, t_end)) +
  scale_x_continuous("Time") +
  geom_vline(xintercept = tq2, linetype = "dotted") +
  annotate("text", x = 1.06*tq2, y = 3, label = "t[iq2]", parse = TRUE, hjust = 1.1, size = 6) +
  geom_segment(
    aes(x = 1.03*tq2, y = 2.9, xend =  1.01*tq2, yend = 2.7),
    arrow = arrow(length = unit(0.3, "cm")),
  ) +
  geom_vline(xintercept = tq1, linetype = "dotted") +
  annotate("text", x = 0.98*tq1, y = 3, label = "t[iq1]", parse = TRUE, hjust = 1.1, size = 6) +
  geom_segment(
    aes(x = 0.97*tq1, y = 2.8, xend =  0.99*tq1, yend = 2.6),
    arrow = arrow(length = unit(0.3, "cm")),
  ) +
  geom_vline(xintercept = (Ti_est), linetype = 2) +
  annotate("text", x = 0.985*(Ti_est), y = 3, label = expression(hat(T)[JM]), parse = TRUE, hjust = 1.1, size = 6) +
  geom_segment(
    aes(x = 0.98*(Ti_est), y = 2.8, xend =  0.99*(Ti_est), yend = 2.6),
    arrow = arrow(length = unit(0.3, "cm")),
  ) +
  geom_vline(xintercept = subTidata$Ti_true, linetype = 2) +
  annotate("text", x = subTidata$Ti_true * 1.05, y = 3, label = expression(T[i]), parse = TRUE, hjust = 1.1, size = 6) +
  geom_segment(
    aes(x = 1.03*subTidata$Ti_true, y = 2.9, xend =  1.01*subTidata$Ti_true, yend = 2.7),
    arrow = arrow(length = unit(0.3, "cm")),
  ) +
  geom_hline(yintercept = censor_value, linetype = "dotted") +
  annotate("text", y = 1.1 * censor_value, x = 8,
           label = "Detection~Limit", parse = TRUE, hjust = 1.1, size = 6) +
  theme_classic(base_size = 20) +
  scale_linetype_manual(name = "Method", values = c("solid", "dashed")) +
  ylab(bquote("Viral load (in" ~ log[10] ~ "-scale)")) +
  annotate("rect", xmin = Ti_est_lower, xmax = Ti_est_upper, ymin = 0, ymax = max(subdata$y_decay_halfDL), 
           alpha = 0.3) + ylim(0, max(subdata$y_decay_halfDL)) + 
  geom_curve(aes(x = 2.2, y = 1.596, xend = Ti_est, yend = 1.13), curvature = 0.1)
ggsave(paste0("Plots/estimated change points with CI v2", ".png"), width = 20, height = 6)


  coord_cartesian(xlim = c(4, 4.1), ylim = c(1.5, 2.2))
  # geom_function(fun = function(x) (fixed_est$eta1 + disp_est$Lambda.eta1 * randef$tau1i) +
  #                 (fixed_est$eta3) * exp(- (fixed_est$eta2 + disp_est$Lambda.eta2 * randef$tau2i) * x),
  #               xlim = c(0, Ti), aes(color = "Fitted Viral Decay Curve")) +
  # geom_function(fun = function(x) (fixed_est$beta1 + disp_est$Lambda.beta1 * randef$b1i) *
  #                 (x - Ti) / ((x - Ti) + exp(fixed_est$beta2 -
  #                                              (fixed_est$beta3 + disp_est$Lambda.beta3 * randef$b3i) * (x - Ti))) +
  #                 (fixed_est$beta4),
  #               xlim = c(Ti, max(subdata$time_decay)), aes(color = "Fitted Viral Rebound Curve")) +
  # scale_color_manual(values = c("Fitted Viral Decay Curve" = "blue",
  #                               "Fitted Viral Rebound Curve" = "green"))

ggplot(subdata, aes(x = time_decay, y = y_decay)) +
  geom_point(aes(fill = factor(censor_decay)), size = 2, shape = 21, stroke = 0) +
  geom_line() +
  scale_x_continuous("Time") +
  geom_vline(xintercept = tq2, linetype = 2) +
  annotate("text", x = 1.08*tq2, y = 3, label = "t[iq2]", parse = TRUE, hjust = 1.1) +
  geom_vline(xintercept = tq1, linetype = 2) +
  annotate("text", x = tq1, y = 3, label = "t[iq1]", parse = TRUE, hjust = 1.1) +
  geom_vline(xintercept = Ti_poly_table[1, 2], linetype = 2) +
  annotate("text", x = 1.02* Ti_poly_table[1, 2], y = 2, label = expression(hat(T)[naive]), parse = TRUE, hjust = 1.1) +
  geom_vline(xintercept = (Ti0_est + trt_stop), linetype = 2) +
  annotate("text", x = 1.09*(Ti0_est + trt_stop), y = 2, label = expression(hat(T)[JM]), parse = TRUE, hjust = 1.1) +
  geom_vline(xintercept = subTidata$Ti_true, linetype = 2) +
  annotate("text", x = 0.9996*subTidata$Ti_true, y = 2.5, label = expression(T[i]), parse = TRUE, hjust = 1.1) +
  geom_hline(yintercept = censor_value, linetype = 2) +
  annotate("text", y = 1.03 * censor_value, x = 7,
           label = "Detection~Limit", parse = TRUE, hjust = 1.1) +
  scale_y_continuous(bquote("Viral load (in" ~ log[10] ~ "-scale)")) +
  scale_fill_manual(values = c("black", "red"), labels = c("Observed data", "Censored data")) +
  labs(color = "Curve Type", fill = "Data Type") +
  theme_classic() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"))

NC
table_par <- readRDS("SimResults IV/table_par.RDS")
table_disp <- readRDS("SimResults IV/table_disp.RDS")
table_parsd <- readRDS("SimResults IV/table_parsd.RDS")
table_par_naive <- readRDS("SimResults IV/table_par_naive.RDS")
table_disp_naive <- readRDS("SimResults IV/table_disp_naive.RDS")
table_parsd_naive <- readRDS("SimResults IV/table_parsd_naive.RDS")
# boot_sd <- readRDS("SimResults IV/boot_sd.RDS")
table_Ti_list <- readRDS("SimResults IV/table_Ti_list.RDS")
table_TiCI <- readRDS("SimResults IV/table_TiCI.RDS")

table_par2 <- table_par[c(1:100), ]
table_disp2 <- table_disp[c(1:100), ]
table_parsd2 <- table_parsd[c(1:100), ]
# boot_sd2 <- boot_sd[c(1:21), ]
table_TiCI2 <- table_TiCI[c(1:100), ]
apply(table_TiCI2, 2, mean)

simresult <- simtable(table_par2, table_parsd2, par_names, true_value, table_disp2, disppar_names, disp_value, boot_sd2 = NULL)
saveRDS(simresult, "SimResults IV/result_table.RDS")

table_par_naive2 <- table_par_naive[c(1:100), ]
table_disp_naive2 <- table_disp_naive[c(1:100), ]
table_parsd_naive2 <- table_parsd_naive[c(1:100), ]
# boot_sd_naive2 <- boot_sd_naive[c(1:26), ]
simresult_naive <- simtable(table_par_naive2, table_parsd_naive2, par_names_naive, true_value_naive, table_disp_naive2, disppar_names_naive, disp_value_naive, boot_sd2 = NULL)
saveRDS(simresult_naive, "SimResults IV/result_table_naive.RDS")

simresult <- readRDS("SimResults IV/result_table.RDS")
simresult$Parameter <- c("$\\eta_1$", "$\\eta_2$", "$\\eta_3$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$", "$\\alpha_0$", "$\\alpha_1$", "$\\alpha_2$", "$\\omega_1$", "$\\omega_2$", "$\\Sigma_{2,4}$", "$\\Lambda_{\\tau11}$", "$\\Lambda_{\\tau22}$", "$\\Lambda_{b11}$","$\\Lambda_{b33}$","$\\Lambda_{a11}$", "$\\sigma_1$", "$\\sigma_2$", "$\\sigma_3$", "$\\sigma_4$")
latex_table<-xtable(simresult, type = "latex",align=c("ccccccccc"))
digits(latex_table) <- c(0,2,2,3,3,3,3,3,2)
print(latex_table, file = "SimResults IV/latex_table.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,length(true_value), length(simresult$Parameter)))

simresult_naive <- readRDS("SimResults IV/result_table_naive.RDS")
simresult_naive$Parameter <- c("$\\eta_1$", "$\\eta_2$", "$\\eta_3$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$", "$\\omega_1$", "$\\omega_2$", "$\\Sigma_{2,4}$", "$\\Lambda_{\\tau11}$", "$\\Lambda_{\\tau22}$", "$\\Lambda_{b11}$","$\\Lambda_{b33}$", "$\\sigma_1$", "$\\sigma_2$", "$\\sigma_4$")
latex_table<-xtable(simresult_naive, type = "latex",align=c("ccccccccc"))
digits(latex_table) <- c(0,2,2,3,3,3,3,3,2)
print(latex_table, file = "SimResults IV/latex_table_naive.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,length(true_value_naive), length(simresult$Parameter_naive)))


