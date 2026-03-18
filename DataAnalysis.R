library(nlme)
library(tidyverse)
library(optimx)
library(mvtnorm)
library(Matrix)
library(xtable)

setwd(here::here())
setwd("R")
file.sources = list.files(pattern = "*.R")
sapply(file.sources, source, .GlobalEnv)

setwd(here::here())
data_decay <- read.csv("Data/cleaned_data_decay.csv")
data_rebound <- read.csv("Data/cleaned_data_rebound.csv")
data_cd4 <- read.csv("Data/cleaned_data_cd4.csv")
data_trans <- read.csv("Data/cleaned_data_trans.csv")

# Summary
n <- length(unique(data_decay$PATIENT))
followup <- data_decay %>% 
  group_by(PATIENT) %>% 
  reframe(n = n(), t_max = max(time_decay))
followup_CD4 <- data_cd4 %>% 
  group_by(PATIENT) %>% 
  reframe(n = n())

print(paste0("The dataset includes ", n,
             " patients who were followed over time. Viral loads and CD4 counts were repeatedly measured, with follow-up periods ranging from ",
             round(min(followup$t_max),1) , " to ",  round(max(followup$t_max),1), " years. The number of repeated viral load measurements varied across patients, ranging from ", 
             min(followup$n) , " to ",  max(followup$n), " with an average of ",  round(mean(followup$n),1), " measurements per patient. Approximately ",
             round(mean(data_decay$censor_decay)*100,1), "% of viral load measurements are left-censored due to a detection limit of 1.6 on the log10 scale. For CD4 counts, the number of repeated measurements ranged from ",
             min(followup_CD4$n) , " to ",  max(followup_CD4$n), " with an average of ",  mean(followup_CD4$n), " measurements per patient."))

# Viral load trajectories for all patients
summary_by_patient <- data_decay %>%
  group_by(PATIENT) %>%
  reframe(
    last_decay = max(time_decay[treatment == 1], na.rm = TRUE)
  )
data_plot <- data_decay %>%
  left_join(summary_by_patient, by = "PATIENT") %>%
  mutate(treatment = ifelse(time_decay == last_decay, 0, treatment))
ggplot(data_plot, aes(x=time_decay, y=y_decay)) + 
  geom_point(aes(fill=factor(censor_decay)),size=2, shape=21, stroke=0) +
  geom_line(aes(group=PATIENT,color=factor(treatment))) +
  scale_x_continuous("Time (in years)") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
  scale_color_manual(values=c("steelblue","black"),labels = c("Following ART interruption", "During ART"))+
  labs(color = "ART status", fill = "Data type")+ggtitle("Plot for all observations")+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
# ggsave(paste0("Plots/Viral load.png"), width = 10, height = 6)
ggplot(data_plot, aes(x=time_decay, y=y_decay_halfDL)) + 
  geom_point(aes(shape=factor(censor_decay)),size=2) +
  geom_line(aes(group=PATIENT))+
  geom_hline(yintercept=log10(40), linetype="dotted") + 
  annotate("text", x = 5, y=1.8, label = "Detection Limit", size = 5) +
  scale_x_continuous("Time (in years)") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_shape_manual(name = "Censored", values = c(16, 1), labels=c('No', 'Yes'))+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
# ggsave(paste0("Plots/Viral load paper.png"), width = 10, height = 5)
ann_data <- data_frame(
  PATIENT = c(1, 8, 11),
  label = c("Stop Time 1", "Stop Time 2", "Stop Time 3"),
  treatment_stop = c(1.643836, 3.767123, 1.586301),
  x_pos = c(2, 3.4, 1.2), 
  y_pos = c(7, 7, 7)        
)

ann_data2 <- data_plot %>%
  filter(PATIENT %in% c(1, 8, 11)) %>%
  group_by(PATIENT) %>%
  filter(time_decay == max(time_decay)) %>%
  mutate(
    # Mapping the PATIENT IDs to simple 1, 2, 3 labels
    label = case_when(
      PATIENT == 1  ~ "Patient 1",
      PATIENT == 8  ~ "Patient 2",
      PATIENT == 11 ~ "Patient 3"
    ),
    # Nudging the label slightly to the right of the last point
    x_pos = time_decay + 0.1,
    y_pos = y_decay_halfDL
  )

ggplot(data_plot %>% filter(PATIENT %in% c(1, 11, 8)), aes(x=time_decay, y=y_decay_halfDL)) + 
  geom_point(aes(shape=factor(censor_decay)),size=2) +
  geom_line(aes(group=PATIENT))+
  geom_vline(aes(xintercept = treatment_stop), 
             linetype = "dashed", color = "darkgrey") +
  geom_text(data = ann_data, aes(x = x_pos, y = y_pos, label = label), size = 4) +
  geom_curve(data = ann_data, 
             aes(x = x_pos + 0.1, y = y_pos - 0.2, xend = treatment_stop, yend = y_pos - 0.5),
             arrow = arrow(length = unit(0.2, "cm")), 
             curvature = -0.2) +
  
  # Patient labels placed near the end of the trajectories
  geom_text(data = ann_data2, aes(x = x_pos - 0.5, y = y_pos + 0.3, label = label), 
            hjust = 0, size = 4, fontface = "bold") +
  
  geom_hline(yintercept=log10(40), linetype="dotted") + 
  annotate("text", x = 4.5, y=1.8, label = "Detection Limit", size = 5) +
  scale_x_continuous("Time (in years)") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_shape_manual(name = "Censored", values = c(16, 1), labels=c('No', 'Yes'))+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
# ggsave(paste0("Plots/Viral load paper 5 patients.png"), width = 10, height = 5)

# Determine starting values for decay model
data_decay.long <- data_decay %>% filter(time_decay < treatment_stop)
# data_decay.long <- data_decay %>% filter(time_decay <= (first_reboundobs - last_decayobs)*7/8)
data_decay.long <- groupedData(y_decay_halfDL ~ time_decay | PATIENT , data = data_decay.long)

set.seed(2025)
error_mess <- "try-error"
while ("try-error" %in% error_mess){
  start0 <- rnorm(3, c(1.3, 11, 3.8), c(0.5, 1, 0.5))
  model_decay <- try(nlme(y_decay_halfDL ~ eta1 + eta3 * exp(- eta2 * time_decay),
                            data = data_decay.long,
                            fixed = eta1 + eta2 + eta3 ~ 1,
                            random = eta1 + eta2 ~ 1,
                            start = c(eta1 = 1.3, eta2 = 11, eta3 = 3.9),
                            control = nlmeControl(msVerbose = T)))
  error_mess <- attr(model_decay, "class")
}
summary(model_decay)

# model_decay <- nls(y_decay_halfDL ~ eta1 + eta3 * exp(- eta2 * time_decay),
#                     data = data_decay.long,
#                     start = c(eta1 = 1.3, eta2 = 3, eta3 = 3.9))
# summary(model_decay)

randef_decay <- tibble::rownames_to_column(coef(model_decay), "PATIENT")
data_decay_randef <- merge(data_decay.long, randef_decay, by = "PATIENT")
fixed_decay <- summary(model_decay)$coefficients$fixed
Lambda_decay <- as.numeric(VarCorr(model_decay)[, 2])

# uniqueID_decay <- droplevels(unique(data_decay_randef$PATIENT))
# for(i in uniqueID_decay){
#   subdata <- data_decay_randef %>%
#     filter(PATIENT == i)
#   max.time <- max(subdata$time_decay)
#   subrandef <- randef_decay %>% filter(PATIENT == i)
#   eta1 <- subrandef$eta1
#   eta2 <- subrandef$eta2
#   eta3 <- subrandef$eta3
#   subdata$censor_decay <- factor(subdata$censor_decay, levels = c(0, 1), labels = c("No", "Yes"))
#   ggplot(subdata, aes(time_decay, y_decay_halfDL)) +
#     geom_point(aes(shape = censor_decay), size = 3) +
#     scale_shape_manual(name = "Censored", values = c(16, 1)) +
#     xlim(0, 1.1 * max.time) +
#     geom_function(fun = function(x) eta1 + eta3 * exp(- eta2 * x)) +
#     annotate("text", x = max.time, y = 1.8, label = "detection limit", size = 4) +
#     geom_hline(yintercept = log10(40), linetype = "dotted") +
#     xlab("Time in years") +
#     scale_y_continuous(bquote("Viral load (in" ~ log[10] ~ "-scale)")) +
#     theme_classic() +
#     theme(text = element_text(size = 14))
# 
#   ggsave(paste0("Decay_plot/ID", i, ".png"), width = 10, height = 6)
#   cat("i=", i, '\n')
# }
tau2i <- tibble::rownames_to_column(ranef(model_decay), "PATIENT")[, c(1,3)]
tau2i$eta2 <- scale(tau2i$eta2, center = TRUE, scale = TRUE)

# Determine starting values for rebound model
data_rebound.long <- data_rebound %>% 
  # mutate(cutoff = (first_reboundobs - last_decayobs)*7/8) %>% 
  mutate(cutoff = treatment_stop) %>% 
  filter(time_rebound >= cutoff) %>%
  mutate(time_rebound = time_rebound - cutoff)
data_rebound.long <- groupedData(y_rebound ~ time_rebound | PATIENT , data = data_rebound.long)

set.seed(1)
error_mess <- "try-error"
while ("try-error" %in% error_mess){
  start0 <- rnorm(5, c(2.3, 3.2, 1, 1.4, 0.1), c(0.1, 0.1, 0.1, 0.1, 0.1))
  model_rebound <- try(nlme(y_rebound ~ beta1 * time_rebound / (time_rebound + exp(beta2 - beta3 * time_rebound)) + beta4/(1 + exp(beta5 * time_rebound)),
                            data = data_rebound.long,
                            fixed = list(beta1 + beta2 + beta3 + beta4 + beta5 ~ 1),
                            random = beta1 + beta3 ~ 1,
                            start = start0,
                            control = nlmeControl(msVerbose = T)))
  # start0 <- rnorm(4, c(2.3, 3.2, 1, 1.4), c(0.1, 0.1, 0.1, 0.1))
  # model_rebound <- try(nlme(y_rebound ~ beta1 * time_rebound / (time_rebound + exp(beta2 - beta3 * time_rebound)) + beta4,
  #                           data = data_rebound.long,
  #                           fixed = list(beta1 + beta2 + beta3 + beta4 ~ 1),
  #                           random = beta1 + beta3 + beta4 ~ 1,
  #                           start = start0,
  #                           control = nlmeControl(msVerbose = T)))
  error_mess <- attr(model_rebound, "class")
}
summary(model_rebound)

randef <- tibble::rownames_to_column(coef(model_rebound), "PATIENT")
data_rebound_randef <- merge(data_rebound.long, randef, by = "PATIENT")
fixed_rebound <- summary(model_rebound)$coefficients$fixed
Lambda_rebound <- as.numeric(VarCorr(model_rebound)[, 2])

uniqueID_rebound <- unique(data_rebound_randef$PATIENT)
# for(i in uniqueID_rebound){
#   subdata <- data_rebound_randef %>%
#     filter(PATIENT == i)
#   max.time <- max(subdata$time_rebound)
#   subdata$censor_rebound <- factor(subdata$censor_rebound, levels = c(0, 1), labels = c("No", "Yes"))
#   subrandef <- randef %>% filter(PATIENT == i)
#   beta1 <- subrandef$beta1
#   beta2 <- subrandef$beta2
#   beta3 <- subrandef$beta3
#   beta4 <- subrandef$beta4
#   beta5 <- subrandef$beta5
#   ggplot(subdata, aes(time_rebound, y_rebound_halfDL)) +
#     geom_point(aes(shape = censor_rebound), size = 3) +
#     scale_shape_manual(name = "Censored", values = c(16, 1)) +
#     xlim(0, 1.1 * max.time) +
#     geom_function(fun = function(x)  beta1 * x / (x + exp(beta2 - beta3 * x)) + beta4/(1 + exp(beta5 * x))) +
#     # geom_function(fun = function(x)  beta1 * x / (x + exp(beta2 - beta3 * x)) + beta4) +
#     annotate("text", x = max.time, y = 1.8, label = "detection limit", size = 4) +
#     geom_hline(yintercept = log10(40), linetype = "dotted") +
#     xlab("Time in years") +
#     scale_y_continuous(bquote("Viral load (in" ~ log[10] ~ "-scale)")) +
#     theme_classic() +
#     theme(text = element_text(size = 14))
# 
#   ggsave(paste0("Rebound_plot/ID", i, ".png"), width = 10, height = 6)
#   cat("i=", i, '\n')
# }

# Determine starting values for CD4 model
data_cd4.long <- groupedData(y_logCD4 ~ time_CD4 | PATIENT , data = data_cd4)
model_cd4 <- lme(y_logCD4 ~ time_CD4 + I(time_CD4^2),
                 data = data_cd4.long,
                 random = ~ 1 | PATIENT)
summary(model_cd4)

data_cd4_randef <- tibble::rownames_to_column(coef(model_cd4), "PATIENT")
colnames(data_cd4_randef)[2:4] <- c("alpha0", "alpha1", "alpha2")
fixed_cd4 <- summary(model_cd4)$coefficients$fixed
Lambda_cd4 <- as.numeric(VarCorr(model_cd4)[, 2])

uniqueID_cd4 <- unique(data_cd4_randef$PATIENT)

# for(i in uniqueID_decay){
#   subdata <- data_cd4 %>%
#     filter(PATIENT == i)
#   subdata.randef <- data_cd4_randef %>%
#     filter(PATIENT == i)
#   max.time <- max(subdata$time_CD4)
#   ggplot(subdata, aes(time_CD4, y_logCD4)) +
#     geom_point(size = 3) +
#     xlim(0, 1.1 * max.time) +
#     geom_function(fun = function(x) subdata.randef$alpha0 + subdata.randef$alpha1 * x + subdata.randef$alpha2 * x ^ 2) +
#     xlab("Time in years") +
#     scale_y_continuous(bquote("CD4 cell count (in" ~ log[10] ~ "-scale cells/" ~ mm^3 ~ ")")) +
#     theme_classic() +
#     theme(text = element_text(size = 14))
# 
#   ggsave(paste0("CD4_plot/ID", i, ".png"), width = 10, height = 6)
#   cat("i=", i, '\n')
# }

maxCD4 <- merge(data_cd4.long, data_cd4_randef) %>% 
  mutate(CD4est = alpha0 + alpha1 * time_CD4 + alpha2 * time_CD4 ^ 2) %>% 
  group_by(PATIENT) %>% 
  filter(CD4est == max(CD4est)) %>% 
  ungroup()
maxCD4 <- merge(maxCD4, tau2i) %>% dplyr::select(PATIENT, CD4est, eta2)
colnames(maxCD4)[2:3] <- c("maxCD4", "tau2i")

data_trans <- merge(data_trans, maxCD4, all.x = TRUE) %>% 
  mutate(y_guess = (first_reboundobs + treatment_stop)/2)
model_trans <- lm(log(y_guess) ~ tau2i + maxCD4 - 1, data = data_trans)
summary(model_trans)
fixed_trans <- summary(model_trans)$coefficients[, 1]
sigma_trans <- summary(model_trans)$sigma

# Form models
Object_decay <- list(
  modeltype = "nlme",
  response = "y_decay",
  reg_equation = "(eta1 + Lambda.eta1 * tau1i) + eta3 * exp(- (eta2 + Lambda.eta2 * tau2i) * time_decay)",
  distribution = 'normal',
  random_effects = c("tau1i", "tau2i"),
  fixed_param = list(names = c("eta1", "eta2", "eta3"),
                     start_values = c(fixed_decay[1], fixed_decay[2], fixed_decay[3]),
                     # start_values = c(0.92, 10.8, 4),
                     lower = c(-Inf, -Inf, -Inf),
                     upper = c( Inf, Inf, Inf)),
  disp_param = list(names = c("Lambda.eta1", "Lambda.eta2"),
                    start_values = c(Lambda_decay[1], Lambda_decay[2]),
                    lower = c(0, 0),
                    upper = c(Inf, Inf)),
  sigma = Lambda_decay[3],
  left_censoring = list(indicator = "censor_decay",
                        limit_value = log10(40),
                        method = "tobit"))

Object_rebound <- list(
  modeltype = "nlme",
  response = "y_rebound",
  reg_equation = "(beta1 + Lambda.beta1 * b1i) * time_rebound / (time_rebound + exp(beta2  - (beta3+ Lambda.beta3 * b3i) * time_rebound)) + (beta4)/(1 + exp((beta5)* time_rebound))",
  # reg_equation = "(beta1 + Lambda.beta1 * b1i) * time_rebound / (time_rebound + exp(beta2 - (beta3 + Lambda.beta3 * b3i) * time_rebound)) + beta4",
  distribution = 'normal',
  random_effects = c("b1i", "b3i"),
  fixed_param = list(names = c("beta1","beta2", "beta3", "beta4", "beta5"),
                     # start_values = c(fixed_rebound[1], fixed_rebound[2], fixed_rebound[3], fixed_rebound[4]),
                     start_values = c(fixed_rebound[1], fixed_rebound[2], fixed_rebound[3], fixed_rebound[4], fixed_rebound[5]),
                     # start_values = c(3.5, -0.4, 5.6, 0.9, -1), 
                     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf),
                     upper = c( Inf, Inf, Inf, Inf, Inf)),
  disp_param = list(names = c("Lambda.beta1", "Lambda.beta3"),
                    start_values = c(Lambda_rebound[1], Lambda_rebound[2]),
                    lower = c(0, 0),
                    upper = c(Inf, Inf)),
  sigma = Lambda_rebound[3],
  left_censoring = list(indicator = "censor_rebound",
                        limit_value = log10(40),
                        method = "tobit"))

Object_CD4 <- list(
  modeltype = "nlme",
  response = "y_logCD4",
  # reg_equation = "(alpha0 + Lambda.alpha0 * a0i) + (alpha1 + Lambda.alpha1 * a1i) * time_CD4 + (alpha2 + Lambda.alpha2 * a2i) * time_CD4 ^ 2",
  reg_equation = "(alpha0 + Lambda.alpha0 * a0i) + (alpha1) * time_CD4 + (alpha2) * time_CD4 ^ 2",
  distribution = 'normal',
  random_effects = c('a0i'),
  fixed_param = list(names = paste0("alpha", 0:2),
                     start_values = c(fixed_cd4[1], fixed_cd4[2], fixed_cd4[3]),
                     lower = c(-Inf, -Inf, -Inf),
                     upper = c( Inf, Inf, Inf)),
  # disp_param = list(names = c("Lambda.alpha0", "Lambda.alpha1", "Lambda.alpha2"),
  #                   start_values = c(Lambda_cd4[1], Lambda_cd4[2], Lambda_cd4[3]),
  #                   lower = c(0, 0, 0),
  #                   upper = c(Inf, Inf, Inf)),
  disp_param = list(names = c("Lambda.alpha0"),
                    start_values = c(Lambda_cd4[1]),
                    lower = c(0),
                    upper = c(Inf)),
  sigma = Lambda_cd4[2])

Object_trans <- list(
  modeltype = "survival",
  response = "time",
  event = "trans",
  reg_equation = "omega1 * tau2i + omega2 * maxCD4",
  distribution = 'log-normal',
  fixed_param = list(names = c("omega1", "omega2"),
                     # start_values = c(fixed_trans[1], fixed_trans[2]),
                     start_values = c(0.36, -0.37), 
                     lower = c(-Inf, -Inf),
                     upper = c(Inf, Inf)),
  sigma = sigma_trans,
  truncated = list(lower = "treatment_stop",
                   upper = "first_reboundobs"))


memList <- list(Object_decay, Object_rebound, Object_CD4, Object_trans)
dataList <- list(data_decay, data_rebound, data_cd4, data_trans)
dataList <- lapply(dataList, function(df) {
  df$last_decayobs <- 0
  df
})


i = 1
error_mess = "try-error"
while ("try-error" %in% error_mess){
  set.seed(i)
  print(glue::glue("i = ",i))
  md0 <- try(fit_multi_mems(memList = memList,
                            dataList = dataList,
                            subject_id = "PATIENT",
                            randeff_info = list(distribution = "normal", degree = 0),
                            cov_method = "spherical",
                            loglike_tol = 5e-2, par_tol = 5e-2,
                            Silent = T, iterMax = 30, REML = FALSE, naive = FALSE, adjust = FALSE), silent = T)
  error_mess <- attr(md0, "class")
  i = i + 1
}

set.seed(1)
md0 <- try(fit_multi_mems(memList = memList,
                          dataList = dataList,
                          subject_id = "PATIENT",
                          randeff_info = list(distribution = "normal", degree = 0),
                          cov_method = "spherical",
                          loglike_tol = 5e-2, par_tol = 5e-2,
                          Silent = T, iterMax = 30, REML = FALSE, naive = FALSE, adjust = FALSE), silent = T)

fixed_est <- as.list(md0$fixed_est)
disp_est <- as.list(md0$disp_est)
par <- c("$\\eta_1$", "$\\eta_2$", "$\\eta_3$", 
         "$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$","$\\beta_5$",
         "$\\alpha_0$", "$\\alpha_1$", "$\\alpha_2$", 
         "$\\omega_1$", "$\\omega_2$")
disp <- c("$\\Lambda_{\\tau11}$", "$\\Lambda_{\\tau22}$", "$\\Lambda_{b11}$", "$\\Lambda_{b33}$", 
          "$\\Lambda_{a00}$", 
          "$\\sigma_1$", "$\\sigma_2$", "$\\sigma_3$", "$\\sigma_4$")
npar <- length(par)
ndisp <- length(disp)
result_table <- data.frame(matrix(NA, nrow = npar + ndisp, ncol = 5))
colnames(result_table) <- c("Parameter", "Estimate", "SE", "$z$-value", "$p$-value")
result_table$Parameter <- c(par, disp)
result_table$Estimate <- c(fixed_est, disp_est)
result_table$SE <- c(md0$fixed_sd, rep(NA, ndisp))
result_table$`$z$-value` <- as.numeric(result_table$Estimate) / as.numeric(result_table$SE)
result_table$`$p$-value` <- 2 * pnorm(abs(result_table$`$z$-value`), lower.tail = FALSE)
latex_table <- xtable(result_table, type = "latex",align=c("cccccc"))
digits(latex_table) <- c(0,3,3,3,3,3)
print(latex_table, file = "DataAnalysisResultsv2.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,npar+ndisp))

cor_table <- solve(md0$invSIGMA)
round(cor_table,3)

# Fitted plot
uniqueID <- unique(data_decay$PATIENT)

for(i in uniqueID){
  subdat <- data_decay %>% filter(PATIENT == i)
  ranef <- dplyr::filter(md0$Bi, PATIENT == i)
  subdat_cd4 <- merge(dplyr::filter(data_cd4, PATIENT == i), randef) %>% 
    # mutate(CD4est = (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i) + (fixed_est$alpha1 + disp_est$Lambda.alpha1 * ranef$a1i) * time_CD4 + (fixed_est$alpha2 + disp_est$Lambda.alpha1 * ranef$a2i) * time_CD4 ^ 2)
    mutate(CD4est = (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i) + (fixed_est$alpha1) * time_CD4 + (fixed_est$alpha2) * time_CD4 ^ 2)
  subdat_cd4$y_cd4_recale <- rescale(subdat_cd4$y_logCD4, 
                                     from = c(min(subdat_cd4$y_logCD4), max(subdat_cd4$y_logCD4)), 
                                     to = c(min(subdat$y_decay_halfDL), max(subdat$y_decay_halfDL)))
  alpha <- (max(subdat$y_decay_halfDL) - min(subdat$y_decay_halfDL))/(max(subdat_cd4$y_logCD4) - min(subdat_cd4$y_logCD4))
  t_end <- max(subdat$time_decay)
  treatment_stop <- unique(subdat$treatment_stop)
  Ti_est <- md0$Ti %>% filter(PATIENT == i) %>% pull(Ti0_sim) + treatment_stop
  subdat_decay <- subdat %>% filter(time_decay <= Ti_est)
  subdat_rebound <- subdat %>% filter(time_decay > Ti_est) %>% 
    mutate(time_rebound = time_decay - Ti_est)
  upper_bound <- unique(subdat_decay$first_reboundobs)
  top <- 1.1*max(subdat$y_decay_halfDL)

  p <- ggplot() +
    xlim(0, t_end) +
    ylim(0, top) + ggtitle(paste("Paitient ID = ", i)) +
    # change background colors
    geom_rect(aes(ymin=0, ymax=top, xmin=0, xmax=treatment_stop, fill="during ART"), alpha = .8) +
    geom_rect(aes(ymin=0, ymax=top, xmin=treatment_stop, xmax = t_end, fill="after ART stop"), alpha = .8) +
    scale_fill_brewer(palette = 'Pastel2', name = 'Treatment Stage') +
    theme_classic() +
    # plot Ti est
    geom_vline(xintercept = Ti_est, linetype = 2) +
    annotate("text", x = 0.98 * Ti_est, y = 0.8*top, label = expression(hat(T)[JM]), parse = TRUE, hjust = 1.1) +
    # plot data
    geom_point(data = subdat, aes(x = time_decay, y = y_decay_halfDL, shape = factor(censor_decay))) +
    # geom_point(data = subdat_cd4, aes(x = time_CD4, y = y_cd4_recale, colour = "log(CD4)")) +
    # geom_function(fun = function(x) (alpha * (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i - min(subdat_cd4$y_logCD4)) + min(subdat$y_decay_halfDL) + 
    #                                    (fixed_est$alpha1 * alpha) * x + (fixed_est$alpha2 * alpha) * x ^ 2), 
    #               aes(color = "log(CD4)")) + 
    geom_function(fun = function(x) (fixed_est$eta1 + disp_est$Lambda.eta1 * ranef$tau1i) + fixed_est$eta3 * exp(- (fixed_est$eta2 + disp_est$Lambda.eta2 * ranef$tau2i) * x), 
                   xlim = c(0, Ti_est)) + 
    geom_function(fun = function(x) ((fixed_est$beta1 + disp_est$Lambda.beta1 * ranef$b1i) * (x - Ti_est) / ((x - Ti_est) + exp(fixed_est$beta2 - (fixed_est$beta3+ disp_est$Lambda.beta3 * ranef$b3i) * (x - Ti_est))) + (fixed_est$beta4 )/(1+exp((fixed_est$beta5)* (x - Ti_est)))),
                  xlim = c(Ti_est, t_end)) +
    # geom_function(fun = function(x) ((fixed_est$beta1 + disp_est$Lambda.beta1 * ranef$b1i) * (x - Ti_est) / ((x - Ti_est) + exp(fixed_est$beta2 - (fixed_est$beta3  + disp_est$Lambda.beta3 * ranef$b3i) * (x - Ti_est))) + fixed_est$beta4 ),
    #               xlim = c(Ti_est, t_end)) +
    geom_hline(yintercept = log10(40), linetype = 'dotted') +
    geom_text(data = subdat, x = 9/10*t_end, y = 1.4, label = "detection limit") +
    geom_vline(xintercept = upper_bound, linetype = 'dotted') +
    geom_text(data = subdat, x = 1.05*upper_bound, y = 3/4 * max(subdat$y_decay_halfDL), label = "Latest\nchange\npoint") +
    scale_shape_manual(name = "Censor", values=c(19, 1), labels = c("No", "Yes")) +
    xlab("Time in days") + ylab(bquote("Viral load (in" ~ log[10]~"-scale)"))
  p
  ggsave(paste0("Ind_plot_fitted_paper_v2/ID", i, " Indfit", ".png"), width = 10, height = 6)
  cat("i=", i, '\n')
}

Ti_table <- data_cd4 %>%
  left_join(md0$Bi, by = "PATIENT") %>%
  left_join(md0$Ti, by = "PATIENT") %>% 
  # mutate(fitted_CD4 = fixed_est$alpha0 + disp_est$Lambda.alpha0 * a0i + (fixed_est$alpha1 + disp_est$Lambda.alpha1 * a1i) * time_CD4 + (fixed_est$alpha2 + disp_est$Lambda.alpha2 * a2i) * time_CD4 ^ 2) %>% 
  mutate(fitted_CD4 = fixed_est$alpha0 + disp_est$Lambda.alpha0 * a0i + (fixed_est$alpha1) * time_CD4 + (fixed_est$alpha2) * time_CD4 ^ 2) %>% 
  group_by(PATIENT) %>% 
  mutate(max_CD4 = max(fitted_CD4)) %>% 
  dplyr::select(-c(time_CD4, y_CD4, fitted_CD4, y_logCD4)) %>% 
  distinct() %>% 
  mutate(Ti_est = exp(fixed_est$omega1 * tau2i + fixed_est$omega2 * max_CD4))

# Naive (polynomial) method
# data_poly <- data_decay %>%
#   arrange(PATIENT, time_decay) %>%
#   group_by(PATIENT) %>%
#   mutate(
#     rebound_count = cumsum(treatment == 0 & censor_decay == 0)
#   ) %>%
#   filter(rebound_count <= 5 | treatment == 1) %>% 
#   mutate(y_poly = ifelse(censor_decay == 1, log10(20), y_decay))
# data_poly <- groupedData(y_decay_halfDL ~ time_decay | PATIENT, data = data_decay)
# model_poly <- lme(y_decay_halfDL ~ time_decay + I(time_decay^2), data = data_poly, 
#                   random = ~ time_decay | PATIENT, method = "ML")
# # Fixed effects
# a_fixed <- fixed.effects(model_poly)["(Intercept)"]
# b_fixed <- fixed.effects(model_poly)["time_decay"]
# c_fixed <- fixed.effects(model_poly)["I(time_decay^2)"]
# # Extract random effects
# ranefs <- ranef(model_poly)  # Data frame with rownames = patient IDs
# aibi <- rownames_to_column(ranefs, "ID")
# Ti_poly <- data.frame(ID = rownames(ranefs), Ti = NA)
# for (id in rownames(ranefs)) {
#   # Patient-specific b_i
#   b_i <- b_fixed + aibi[which(aibi$ID == id), ]$time_decay
#   
#   # 95% CI
#   Ti_poly[Ti_poly$ID == id, "Ti"] <- -b_i / (2 * c_fixed)
# }


# Individual fitted plot
# uniqueID <- unique(data_decay$PATIENT)
# 
# for(i in uniqueID){
#   subdat <- data_decay %>% filter(PATIENT == i)
#   ranef <- dplyr::filter(md0$Bi, PATIENT == i)
#   subdat_cd4 <- merge(dplyr::filter(data_cd4, PATIENT == i), randef) %>% 
#     # mutate(CD4est = (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i) + (fixed_est$alpha1 + disp_est$Lambda.alpha1 * ranef$a1i) * time_CD4 + (fixed_est$alpha2 + disp_est$Lambda.alpha1 * ranef$a2i) * time_CD4 ^ 2)
#     mutate(CD4est = (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i) + (fixed_est$alpha1) * time_CD4 + (fixed_est$alpha2) * time_CD4 ^ 2)
#   subdat_cd4$y_cd4_recale <- rescale(subdat_cd4$y_logCD4, 
#                                      from = c(min(subdat_cd4$y_logCD4), max(subdat_cd4$y_logCD4)), 
#                                      to = c(min(subdat$y_decay_halfDL), max(subdat$y_decay_halfDL)))
#   alpha <- (max(subdat$y_decay_halfDL) - min(subdat$y_decay_halfDL))/(max(subdat_cd4$y_logCD4) - min(subdat_cd4$y_logCD4))
#   t_end <- max(subdat$time_decay)
#   Ti_est <- Ti_table %>% filter(PATIENT == i) %>% pull(Ti0_sim) 
#   subdat_decay <- subdat %>% filter(time_decay <= Ti_est)
#   subdat_rebound <- subdat %>% filter(time_decay > Ti_est) %>% 
#     mutate(time_rebound = time_decay - Ti_est)
#   lower_bound <- unique(subdat_decay$last_decayobs)
#   upper_bound <- unique(subdat_decay$first_reboundobs)
#   treatment_stop <- unique(subdat_decay$treatment_stop)
#   Ti_poly_est <- Ti_poly[Ti_poly$ID == i, "Ti"]
#   top <- 1.1*max(subdat$y_decay_halfDL)
#   ai <- a_fixed + aibi[which(aibi$ID == i), ]$`(Intercept)`
#   bi <- b_fixed + aibi[which(aibi$ID == i), ]$time_decay
#   
#   
#   p <- ggplot() +
#     xlim(0, t_end) +
#     ylim(0, top) + ggtitle(paste("Paitient ID = ", i)) +
#     # change background colors
#     geom_rect(aes(ymin=0, ymax=top, xmin=0, xmax=treatment_stop, fill="during ART"), alpha = .8) +
#     geom_rect(aes(ymin=0, ymax=top, xmin=treatment_stop, xmax = t_end, fill="after ART stop"), alpha = .8) +
#     scale_fill_brewer(palette = 'Pastel2', name = 'Treatment Stage') +
#     theme_classic() +
#     # plot Ti est
#     geom_vline(xintercept = Ti_est, linetype = 2) +
#     annotate("text", x = 0.98 * Ti_est, y = 4, label = expression(hat(T)[JM]), parse = TRUE, hjust = 1.1) +
#     geom_vline(xintercept = Ti_poly_est, linetype = 2) +
#     annotate("text", x = 0.98 * Ti_poly_est, y = 4, label = expression(hat(T)[naive]), parse = TRUE, hjust = 1.1) +
#     # plot data
#     geom_point(data = subdat, aes(x = time_decay, y = y_decay_halfDL, shape = factor(censor_decay))) +
#     # geom_point(data = subdat_cd4, aes(x = time_CD4, y = y_cd4_recale, colour = "log(CD4)")) +
#     # geom_function(fun = function(x) (alpha * (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i - min(subdat_cd4$y_logCD4)) + min(subdat$y_decay_halfDL) + 
#     #                                    (fixed_est$alpha1 * alpha) * x + (fixed_est$alpha2 * alpha) * x ^ 2), 
#     #               aes(color = "log(CD4)")) + 
#     geom_function(fun = function(x) (fixed_est$eta1 + disp_est$Lambda.eta1 * ranef$tau1i) + fixed_est$eta3 * exp(- (fixed_est$eta2 + disp_est$Lambda.eta2 * ranef$tau2i) * x), 
#                   aes(linetype = "JM"), xlim = c(0, Ti_est)) + 
#     geom_function(fun = function(x) ((fixed_est$beta1 + disp_est$Lambda.beta1 * ranef$b1i) * (x - Ti_est) / ((x - Ti_est) + exp(fixed_est$beta2 - (fixed_est$beta3 + disp_est$Lambda.beta3 * ranef$b3i) * (x - Ti_est))) + fixed_est$beta4/(1+exp((fixed_est$beta5 + disp_est$Lambda.beta5 * ranef$b5i)* (x - Ti_est)))),
#                   aes(linetype = "JM"), xlim = c(Ti_est, t_end)) +
#     geom_function(fun = function(x) (ai + bi * x + c_fixed * x^2),
#                   aes(linetype = "Naive"), xlim = c(0, t_end)) +
#     # geom_function(fun = function(x) ((fixed_est$beta1 + disp_est$Lambda.beta1 * ranef$b1i) * (x - Ti_est) / ((x - Ti_est) + exp(fixed_est$beta2 - (fixed_est$beta3 ) * (x - Ti_est))) + fixed_est$beta4  + disp_est$Lambda.beta4 * ranef$b4i ),
#     #               aes(color = "log10(RNA)"), xlim = c(Ti_est, t_end)) +
#     geom_hline(yintercept = log10(40), linetype = 'dotted') +
#     geom_text(data = subdat, x = 9/10*t_end, y = 1.4, label = "detection limit") +
#     geom_vline(xintercept = lower_bound, linetype = 'dotted') +
#     geom_text(data = subdat, x = 1.1*lower_bound, y = 2/3 * max(subdat$y_decay_halfDL), label = "Earlist\nchange\npoint") +
#     geom_vline(xintercept = upper_bound, linetype = 'dotted') +
#     geom_text(data = subdat, x = 1.05*upper_bound, y = 2/3 * max(subdat$y_decay_halfDL), label = "Latest\nchange\npoint") +
#     scale_shape_manual(name = "Censor", values=c(19, 1), labels = c("No", "Yes")) +
#     scale_linetype_manual(name = "Method", values = c("solid", "dashed")) +
#     xlab("Time in days") + ylab(bquote("Viral load (in" ~ log[10]~"-scale)"))
#   p
#   ggsave(paste0("Ind_plot_fitted_paper/ID", i, " Indfit", ".png"), width = 10, height = 6)
#   cat("i=", i, '\n')
# }

i=67
subdat <- data_decay %>% filter(PATIENT == i)
ranef <- dplyr::filter(md0$Bi, PATIENT == i)
subdat_cd4 <- merge(dplyr::filter(data_cd4, PATIENT == i), randef) %>% 
  # mutate(CD4est = (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i) + (fixed_est$alpha1 + disp_est$Lambda.alpha1 * ranef$a1i) * time_CD4 + (fixed_est$alpha2 + disp_est$Lambda.alpha1 * ranef$a2i) * time_CD4 ^ 2)
  mutate(CD4est = (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i) + (fixed_est$alpha1) * time_CD4 + (fixed_est$alpha2) * time_CD4 ^ 2)
subdat_cd4$y_cd4_recale <- rescale(subdat_cd4$y_logCD4, 
                                   from = c(min(subdat_cd4$y_logCD4), max(subdat_cd4$y_logCD4)), 
                                   to = c(min(subdat$y_decay_halfDL), max(subdat$y_decay_halfDL)))
alpha <- (max(subdat$y_decay_halfDL) - min(subdat$y_decay_halfDL))/(max(subdat_cd4$y_logCD4) - min(subdat_cd4$y_logCD4))
treatment_stop <- unique(subdat$treatment_stop)
Ti_est <- md0$Ti %>% filter(PATIENT == i) %>% pull(Ti0_sim) + treatment_stop
subdat_decay <- subdat %>% filter(time_decay <= Ti_est)
subdat_rebound <- subdat %>% filter(time_decay > Ti_est) %>% 
  mutate(time_rebound = time_decay - Ti_est)
upper_bound <- unique(subdat_decay$first_reboundobs)
top <- 1.1*max(subdat$y_decay_halfDL)

subdat <- subdat %>% 
  mutate(time_decay_new = ifelse(time_decay > Ti_est, time_decay + 0.03, time_decay))
t_end <- max(subdat$time_decay_new)
(p <- ggplot() +
    xlim(0, t_end) +
    ylim(0, top) +
    theme_classic() +
    # plot Ti est
    geom_vline(xintercept = Ti_est, linetype = 2) +
    annotate("text", x =1.04 * Ti_est, y = 0.85*top, label = expression(hat(T)[JM]), parse = TRUE, hjust = 1.1, size = 6) +
    geom_segment(
      aes(x = 1.02 * Ti_est, y = 0.8*top, xend =  1.01 * Ti_est, yend = 0.75*top),
      arrow = arrow(length = unit(0.3, "cm")), 
    ) +
    # plot data
    geom_point(data = subdat, aes(x = time_decay_new, y = y_decay_halfDL, shape = factor(censor_decay)), size = 2) +
    # geom_point(data = subdat_cd4, aes(x = time_CD4, y = y_cd4_recale, colour = "log(CD4)")) +
    # geom_function(fun = function(x) (alpha * (fixed_est$alpha0 + disp_est$Lambda.alpha0 * ranef$a0i - min(subdat_cd4$y_logCD4)) + min(subdat$y_decay_halfDL) + 
    #                                    (fixed_est$alpha1 * alpha) * x + (fixed_est$alpha2 * alpha) * x ^ 2), 
    #               aes(color = "log(CD4)")) + 
    geom_function(fun = function(x) (fixed_est$eta1 + disp_est$Lambda.eta1 * ranef$tau1i) + fixed_est$eta3 * exp(- (fixed_est$eta2 + disp_est$Lambda.eta2 * ranef$tau2i) * x), 
                  xlim = c(0, Ti_est)) + 
    geom_function(fun = function(x) ((fixed_est$beta1 + disp_est$Lambda.beta1 * ranef$b1i) * (x - Ti_est-0.03) / ((x - Ti_est-0.03) + exp(fixed_est$beta2 - (fixed_est$beta3+ disp_est$Lambda.beta3 * ranef$b3i) * (x - Ti_est-0.03))) + (fixed_est$beta4 )/(1+exp((fixed_est$beta5)* (x - Ti_est-0.03)))),
                  xlim = c(Ti_est+0.03, t_end)) +
    # geom_function(fun = function(x) ((fixed_est$beta1 + disp_est$Lambda.beta1 * ranef$b1i) * (x - Ti_est) / ((x - Ti_est) + exp(fixed_est$beta2 - (fixed_est$beta3  + disp_est$Lambda.beta3 * ranef$b3i) * (x - Ti_est))) + fixed_est$beta4 ),
    #               xlim = c(Ti_est, t_end)) +
    geom_hline(yintercept = log10(40), linetype = 'dotted') +
    geom_text(data = subdat, x = 9/10*t_end, y = 1.4, label = "Detection limit", size = 6) +
    geom_vline(xintercept = upper_bound, linetype = 'dotted') +
    geom_text(data = subdat, x = 1.08*upper_bound, y = 0.63 * max(subdat$y_decay_halfDL), label = "Latest\nchange\npoint", size = 6) +
    geom_segment(
      aes(x = 1.05*upper_bound, y = 0.48 * max(subdat$y_decay_halfDL), xend =  1.01*upper_bound, yend = 0.4 * max(subdat$y_decay_halfDL)),
      arrow = arrow(),
    ) +
    geom_vline(xintercept = treatment_stop, linetype = 'dotted') +
    geom_text(data = subdat, x = 0.93*treatment_stop, y = 0.63 * max(subdat$y_decay_halfDL), label = "Treatment\nstop\ntime", size = 6) +
    geom_segment(
      aes(x = 0.95*treatment_stop, y = 0.48 * max(subdat$y_decay_halfDL), xend =  0.99*treatment_stop, yend = 0.4 * max(subdat$y_decay_halfDL)),
      arrow = arrow(),
    ) +
    scale_shape_manual(name = "Censored", values=c(19, 1), labels = c("No", "Yes")) +
    xlab("Time in years") + ylab(bquote("Viral load (in" ~ log[10]~"-scale)")) +
    theme_classic(base_size = 20)) + geom_curve(aes(x = Ti_est, y = 0.5, xend = Ti_est+0.03, yend = log10(40)+0.05), curvature = 0.08)
p2 <- p + xlim(1, 1.5) + ylim(0, top)
p | p2
ggsave(paste0("Ind_plot_fitted_paper_v2/ID", i, " Indfit v2", ".png"), width = 20, height = 6)
