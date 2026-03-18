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

table_par <- table_parsd <- boot_sd <- data.frame(matrix(nrow = num_sim, ncol = par_length))
table_disp <- data.frame(matrix(nrow = num_sim, ncol = disppar_length))
colnames(table_par) <- colnames(table_parsd) <- colnames(boot_sd) <- c(par_names)
colnames(table_disp) <- disppar_names

table_TiCI <- data.frame(matrix(nrow = num_sim, ncol = 2))
colnames(table_TiCI) <- c("New", "Bootstrap")
table_Ti_list <- list()

table_par <- readRDS("SimResultsBoot/table_par.RDS")
table_disp <- readRDS("SimResultsBoot/table_disp.RDS")
table_parsd <- readRDS("SimResultsBoot/table_parsd.RDS")
table_Ti_list <- readRDS("SimResultsBoot/table_Ti_list.RDS")
table_TiCI <- readRDS("SimResultsBoot/table_TiCI.RDS")
boot_sd <- readRDS("SimResultsBoot/boot_sd.RDS")


i = 1
NC <- 0
while (i <= 2){
  print(glue::glue("i = ",i))
  
  data_sim <- create_simdata(time, num_patient, censor_value, true_value_par, 
                             disp_value_par, newSigma, trt = T)
  
  
  data_decay <- data_sim$decay
  data_rebound <- data_sim$rebound
  data_CD4 <- data_sim$CD4
  data_trans <- data_sim$trans
  
  mean(data_decay$censor)
  
  # Form models
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
  # saveRDS(memList, "memList.RDS")
  # saveRDS(dataList, "dataList.RDS")
  # memList <- readRDS("memList.RDS")
  # dataList <- readRDS("dataList.RDS")
  
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
  } else if(any(is.na(md0$fixed_sd))){
    NC <- NC + 1
  } else{
    randef <- md0$Bi
    fixed_est <- as.list(md0$fixed_est)
    disp_est <- as.list(c(md0$disp_est, Lambda.zero))
    # Ti <- merge(data_CD4, md0$Bi) %>% 
    #   mutate(fitted_CD4 = fixed_est$alpha0 + disp_est$Lambda.alpha0 * a0i + (fixed_est$alpha1) * time_CD4 + fixed_est$alpha2 * time_CD4 ^ 2) %>% 
    #   group_by(ID) %>% 
    #   mutate(max_CD4 = max(fitted_CD4)) %>% 
    #   dplyr::select(-c(time_CD4, CD4, fitted_CD4)) %>% 
    #   distinct() %>% 
    #   mutate(Ti0 = exp(fixed_est$omega1 * tau2i + fixed_est$omega2 * max_CD4))
    Ti_est <- md0$Ti
    Ti_true <- data_trans %>% 
      mutate(upper_bound = first_reboundobs - treatment_stop) %>% 
      dplyr::select(ID, Ti0, Ti_true, upper_bound)
    colnames(Ti_true)[2] <- "Ti0_true"
    TiCI <- merge(Ti_est, Ti_true) %>% 
      mutate(lowerCI = qlnormTrunc(0.025, meanlog = log(Ti0_sim), sdlog = disp_est$time_sigma, min = 0, max = upper_bound)) %>% 
      mutate(upperCI = qlnormTrunc(0.975, meanlog = log(Ti0_sim), sdlog = disp_est$time_sigma, min = 0, max = upper_bound)) %>% 
      mutate(CI = Ti0_true >= lowerCI & Ti0_true <= upperCI)
    table_par[i, ] <- round(md0$fixed_est, 4)
    table_disp[i, ] <- c(solve(md0$invSIGMA)[2, 3], round(md0$disp_est, 4))
    table_parsd[i, ] <- round(as.numeric(md0$fixed_sd), 4)
    Sigma_est <- matrix(0, nrow = 10, ncol = 10)
    Sigma_est[pos.Lambda.nonzero, pos.Lambda.nonzero] <- solve(md0$invSIGMA)
    diag(Sigma_est) <- 1
    boot_result <- bootFun(time, num_patient, censor_value, par = fixed_est, disp = disp_est, Sigma = Sigma_est, B = 50, i)
    boot_sd[i, ] <- round(apply(boot_result$bootEst, 2, sd), 5)
    TiCI <- TiCI %>% 
      mutate(boot_Ti_sd = round(apply(boot_result$bootTi, 1, sd), 5)) %>% 
      mutate(lowerCI_boot = Ti0_sim - 1.96 * boot_Ti_sd) %>% 
      mutate(upperCI_boot = Ti0_sim + 1.96 * boot_Ti_sd) %>% 
      mutate(CI_boot = Ti0_true >= lowerCI_boot & Ti0_true <= upperCI_boot)
    table_Ti_list[[i]] <- TiCI
    table_TiCI[i, 1] <- mean(TiCI$CI)
    table_TiCI[i, 2] <- mean(TiCI$CI_boot)
    
    saveRDS(table_par, "SimResults IV/table_par.RDS")
    saveRDS(table_disp, "SimResults IV/table_disp.RDS")
    saveRDS(table_parsd, "SimResults IV/table_parsd.RDS")
    saveRDS(table_Ti_list, "SimResults IV/table_Ti_list.RDS")
    saveRDS(table_TiCI, "SimResults IV/table_TiCI.RDS")
    saveRDS(boot_sd, "SimResults IV/boot_sd.RDS")
    i <- i + 1
  }
}

randef <- md0$Bi %>% filter(ID == 1)
fixed_est <- as.list(md0$fixed_est)
disp_est <- as.list(md0$disp_est)
subCD4data <- data_CD4[data_CD4$ID %in% c(1),] %>% 
  mutate(fitted_CD4 = fixed_est$alpha0 + (fixed_est$alpha1 + disp_est$Lambda.alpha1 * randef$a1i) * time_CD4 + fixed_est$alpha2 * time_CD4 ^ 2)
max_CD4 <- max(subCD4data$fitted_CD4)
Ti0 = exp(fixed_est$omega1 * randef$tau2i + fixed_est$omega2 * max_CD4)
Ti = 2.17 + Ti0
subTidata <- data_trans[data_trans$ID %in% c(1),] 

ggplot(subdata, aes(x = time_decay, y = y_decay)) +
  geom_point(aes(fill = factor(censor_decay)), size = 2, shape = 21, stroke = 0) +
  geom_line() +
  scale_x_continuous("Time") +
  geom_vline(xintercept = tq2, linetype = 2) +
  annotate("text", x = 1.1 * tq2, y = 3.2, label = "t[iq2]", parse = TRUE, hjust = 1.1) +
  geom_vline(xintercept = treatmentstop, linetype = 2) +
  annotate("text", x = 0.95 * treatmentstop, y = 3.2, label = "t[iq1]", parse = TRUE, hjust = 1.1) +
  geom_vline(xintercept = Ti, linetype = 2) +
  annotate("text", x = 0.95 * Ti, y = 3.2, label = expression(hat(T)[i]), parse = TRUE, hjust = 1.1) +
  geom_vline(xintercept = subTidata$Ti_true, linetype = 2) +
  annotate("text", x = 0.95 * subTidata$Ti_true, y = 3.2, label = expression(T[i]), parse = TRUE, hjust = 1.1) +
  geom_hline(yintercept = censor_value, linetype = 2) +
  annotate("text", y = 1.1 * censor_value, x = 0.95 * max(subdata$time_decay), 
           label = "Detection~Limit", parse = TRUE, hjust = 1.1) +
  scale_y_continuous(bquote("Viral load (in" ~ log[10] ~ "-scale)")) +
  scale_fill_manual(values = c("black", "red"), labels = c("Observed data", "Censored data")) +
  labs(color = "Curve Type", fill = "Data Type") +
  ggtitle("Viral load trajectory of a typical individual") +
  theme_classic() + 
  theme(text = element_text(size = 14), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90")) +
  geom_function(fun = function(x) (fixed_est$eta1 + disp_est$Lambda.eta1 * randef$tau1i) + 
                  (fixed_est$eta3) * exp(- (fixed_est$eta2 + disp_est$Lambda.eta2 * randef$tau2i) * x), 
                xlim = c(0, Ti), aes(color = "Fitted Viral Decay Curve")) +
  geom_function(fun = function(x) (fixed_est$beta1 + disp_est$Lambda.beta1 * randef$b1i) * 
                  (x - Ti) / ((x - Ti) + exp(fixed_est$beta2 - 
                                               (fixed_est$beta3 + disp_est$Lambda.beta3 * randef$b3i) * (x - Ti))) + 
                  (fixed_est$beta4), 
                xlim = c(Ti, max(subdata$time_decay)), aes(color = "Fitted Viral Rebound Curve")) +
  scale_color_manual(values = c("Fitted Viral Decay Curve" = "blue", 
                                "Fitted Viral Rebound Curve" = "green"))


NC
table_par <- readRDS("SimResultsBootAll/table_par.RDS")
table_disp <- readRDS("SimResultsBootAll/table_disp.RDS")
table_parsd <- readRDS("SimResultsBootAll/table_parsd.RDS")
boot_sd <- readRDS("SimResultsBootAll/boot_sd.RDS")
table_Ti_list <- readRDS("SimResultsBootAll/table_Ti_list.RDS")
table_TiCI <- readRDS("SimResultsBootAll/table_TiCI.RDS")

table_par2 <- na.omit(table_par[c(1:100), ])
table_disp2 <- na.omit(table_disp[c(1:100), ])
table_parsd2 <- na.omit(table_parsd[c(1:100), ])
boot_sd2 <- na.omit(boot_sd[c(1:100), ])
table_TiCI2 <- na.omit(table_TiCI[c(1:100), ])
apply(table_TiCI2, 2, mean)

simresult <- simtable(table_par2, table_parsd2, par_names, true_value, table_disp2, disppar_names, disp_value, boot_sd2)
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
