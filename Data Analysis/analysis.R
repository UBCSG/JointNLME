# =============================================================================
# HIV/AIDS Joint Modeling: Data Analysis
# =============================================================================

library(nlme)
library(tidyverse)
library(optimx)
library(mvtnorm)
library(Matrix)
library(xtable)

setwd(here::here())

# Source all helper functions from the R/ directory
file.sources <- list.files("R", pattern = "\\.R$", full.names = TRUE)
sapply(file.sources, source, .GlobalEnv)


# =============================================================================
# 1. LOAD DATA
# =============================================================================

data_decay   <- read.csv(here::here("Data", "cleaned_data_decay.csv"))
data_rebound <- read.csv(here::here("Data", "cleaned_data_rebound.csv"))
data_cd4     <- read.csv(here::here("Data", "cleaned_data_cd4.csv"))
data_trans   <- read.csv(here::here("Data", "cleaned_data_trans.csv"))


# =============================================================================
# 2. DESCRIPTIVE SUMMARY
# =============================================================================

n          <- length(unique(data_decay$PATIENT))
followup   <- data_decay %>% group_by(PATIENT) %>% reframe(n = n(), t_max = max(time_decay))
followup_cd4 <- data_cd4 %>% group_by(PATIENT) %>% reframe(n = n())

cat(sprintf(
  "The dataset includes %d patients. Follow-up ranged from %.1f to %.1f years. ",
  n, min(followup$t_max), max(followup$t_max)
))
cat(sprintf(
  "Viral load measurements per patient: %d to %d (mean %.1f). ",
  min(followup$n), max(followup$n), mean(followup$n)
))
cat(sprintf(
  "Left-censored viral loads: %.1f%%. ",
  mean(data_decay$censor_decay) * 100
))
cat(sprintf(
  "CD4 measurements per patient: %d to %d (mean %.1f).\n",
  min(followup_cd4$n), max(followup_cd4$n), mean(followup_cd4$n)
))


# =============================================================================
# 3. STARTING VALUES — DECAY MODEL
# =============================================================================
# Fit an nlme model on the ART-on phase to obtain good starting values.

data_decay_on <- data_decay %>%
  filter(time_decay < treatment_stop) %>%
  groupedData(y_decay_halfDL ~ time_decay | PATIENT, data = .)

set.seed(2025)
model_decay <- nlme(y_decay_halfDL ~ eta1 + eta3 * exp(-eta2 * time_decay),
                    data    = data_decay_on,
                    fixed   = eta1 + eta2 + eta3 ~ 1,
                    random  = eta1 + eta2 ~ 1,
                    start   = c(eta1 = 1.3, eta2 = 11, eta3 = 3.9),
                    control = nlmeControl(msVerbose = FALSE))
summary(model_decay)

fixed_decay  <- summary(model_decay)$coefficients$fixed
Lambda_decay <- as.numeric(VarCorr(model_decay)[, 2])

# Standardised random effect for eta2 (used as covariate in transition model)
tau2i <- tibble::rownames_to_column(ranef(model_decay), "PATIENT") %>%
  dplyr::select(PATIENT, eta2) %>%
  mutate(eta2 = as.numeric(scale(eta2)))


# =============================================================================
# 4. STARTING VALUES — REBOUND MODEL
# =============================================================================

data_rebound_off <- data_rebound %>%
  filter(time_rebound >= treatment_stop) %>%
  mutate(time_rebound = time_rebound - treatment_stop) %>%
  groupedData(y_rebound ~ time_rebound | PATIENT, data = .)

set.seed(1)
start0 <- rnorm(5, c(2.3, 3.2, 1, 1.4, 0.1), 0.1)
model_rebound <- nlme(y_rebound ~ beta1 * time_rebound / (time_rebound + exp(beta2 - beta3 * time_rebound)) +
                        beta4 / (1 + exp(beta5 * time_rebound)),
                      data    = data_rebound_off,
                      fixed   = list(beta1 + beta2 + beta3 + beta4 + beta5 ~ 1),
                      random  = beta1 + beta3 ~ 1,
                      start   = start0,
                      control = nlmeControl(msVerbose = FALSE))
summary(model_rebound)

fixed_rebound  <- summary(model_rebound)$coefficients$fixed
Lambda_rebound <- as.numeric(VarCorr(model_rebound)[, 2])


# =============================================================================
# 5. STARTING VALUES — CD4 MODEL
# =============================================================================

data_cd4_long <- groupedData(y_logCD4 ~ time_CD4 | PATIENT, data = data_cd4)

model_cd4 <- lme(y_logCD4 ~ time_CD4 + I(time_CD4^2),
                 data   = data_cd4_long,
                 random = ~1 | PATIENT)
summary(model_cd4)

cd4_randef <- tibble::rownames_to_column(coef(model_cd4), "PATIENT") %>%
  dplyr::rename(alpha0 = `(Intercept)`, alpha1 = time_CD4, alpha2 = `I(time_CD4^2)`)

fixed_cd4  <- summary(model_cd4)$coefficients$fixed
Lambda_cd4 <- as.numeric(VarCorr(model_cd4)[, 2])

# Patient-level fitted peak CD4 value
maxCD4 <- data_cd4_long %>%
  left_join(cd4_randef, by = "PATIENT") %>%
  mutate(CD4est = alpha0 + alpha1 * time_CD4 + alpha2 * time_CD4^2) %>%
  group_by(PATIENT) %>%
  slice_max(CD4est, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(tau2i, by = "PATIENT") %>%
  dplyr::select(PATIENT, maxCD4 = CD4est, tau2i = eta2)
maxCD4$PATIENT <- as.integer(maxCD4$PATIENT )

# =============================================================================
# 6. STARTING VALUES — TRANSITION MODEL
# =============================================================================

data_trans <- data_trans %>%
  left_join(maxCD4, by = "PATIENT") %>%
  mutate(y_guess = (first_reboundobs + treatment_stop) / 2)

model_trans <- lm(log(y_guess) ~ tau2i + maxCD4 - 1, data = data_trans)
summary(model_trans)

fixed_trans <- coef(model_trans)
sigma_trans <- summary(model_trans)$sigma


# =============================================================================
# 7. MODEL SPECIFICATION OBJECTS
# =============================================================================

Object_decay <- list(
  modeltype      = "nlme",
  response       = "y_decay",
  reg_equation   = "(eta1 + Lambda.eta1 * tau1i) + eta3 * exp(-(eta2 + Lambda.eta2 * tau2i) * time_decay)",
  distribution   = "normal",
  random_effects = c("tau1i", "tau2i"),
  fixed_param    = list(names        = c("eta1", "eta2", "eta3"),
                        start_values = fixed_decay,
                        lower        = rep(-Inf, 3),
                        upper        = rep( Inf, 3)),
  disp_param     = list(names        = c("Lambda.eta1", "Lambda.eta2"),
                        start_values = Lambda_decay[1:2],
                        lower        = c(0, 0),
                        upper        = c(Inf, Inf)),
  sigma          = Lambda_decay[3],
  left_censoring = list(indicator    = "censor_decay",
                        limit_value  = log10(40),
                        method       = "tobit")
)

Object_rebound <- list(
  modeltype      = "nlme",
  response       = "y_rebound",
  reg_equation   = "(beta1 + Lambda.beta1 * b1i) * time_rebound / (time_rebound + exp(beta2 - (beta3 + Lambda.beta3 * b3i) * time_rebound)) + beta4 / (1 + exp(beta5 * time_rebound))",
  distribution   = "normal",
  random_effects = c("b1i", "b3i"),
  fixed_param    = list(names        = c("beta1", "beta2", "beta3", "beta4", "beta5"),
                        start_values = fixed_rebound,
                        lower        = rep(-Inf, 5),
                        upper        = rep( Inf, 5)),
  disp_param     = list(names        = c("Lambda.beta1", "Lambda.beta3"),
                        start_values = Lambda_rebound[1:2],
                        lower        = c(0, 0),
                        upper        = c(Inf, Inf)),
  sigma          = Lambda_rebound[3],
  left_censoring = list(indicator    = "censor_rebound",
                        limit_value  = log10(40),
                        method       = "tobit")
)

Object_CD4 <- list(
  modeltype      = "nlme",
  response       = "y_logCD4",
  reg_equation   = "(alpha0 + Lambda.alpha0 * a0i) + alpha1 * time_CD4 + alpha2 * time_CD4^2",
  distribution   = "normal",
  random_effects = "a0i",
  fixed_param    = list(names        = paste0("alpha", 0:2),
                        start_values = fixed_cd4,
                        lower        = rep(-Inf, 3),
                        upper        = rep( Inf, 3)),
  disp_param     = list(names        = "Lambda.alpha0",
                        start_values = Lambda_cd4[1],
                        lower        = 0,
                        upper        = Inf),
  sigma          = Lambda_cd4[2]
)

Object_trans <- list(
  modeltype    = "survival",
  response     = "time",
  event        = "trans",
  reg_equation = "omega1 * tau2i + omega2 * maxCD4",
  distribution = "log-normal",
  fixed_param  = list(names        = c("omega1", "omega2"),
                      start_values = c(0.36, -0.37),
                      lower        = rep(-Inf, 2),
                      upper        = rep( Inf, 2)),
  sigma        = sigma_trans,
  truncated    = list(lower = "treatment_stop",
                      upper = "first_reboundobs")
)


# =============================================================================
# 8. FIT JOINT MODEL
# =============================================================================

memList  <- list(Object_decay, Object_rebound, Object_CD4, Object_trans)
dataList <- list(data_decay, data_rebound, data_cd4, data_trans)

# Ensure last_decayobs column exists in all datasets
dataList <- lapply(dataList, function(df) { df$last_decayobs <- 0; df })

set.seed(1)
md0 <- fit_multi_mems(memList      = memList,
                      dataList     = dataList,
                      subject_id   = "PATIENT",
                      randeff_info = list(distribution = "normal", degree = 0),
                      cov_method   = "spherical",
                      loglike_tol  = 5e-2,
                      par_tol      = 5e-2,
                      Silent       = TRUE,
                      iterMax      = 30,
                      REML         = FALSE,
                      naive        = FALSE,
                      adjust       = FALSE)


# =============================================================================
# 9. RESULTS TABLE (LaTeX)
# =============================================================================

fixed_est <- as.list(md0$fixed_est)
disp_est  <- as.list(md0$disp_est)

par_names  <- c("$\\eta_1$", "$\\eta_2$", "$\\eta_3$",
                "$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$", "$\\beta_5$",
                "$\\alpha_0$", "$\\alpha_1$", "$\\alpha_2$",
                "$\\omega_1$", "$\\omega_2$")
disp_names <- c("$\\Lambda_{\\tau11}$", "$\\Lambda_{\\tau22}$",
                "$\\Lambda_{b11}$",    "$\\Lambda_{b33}$",
                "$\\Lambda_{a00}$",
                "$\\sigma_1$", "$\\sigma_2$", "$\\sigma_3$", "$\\sigma_4$")

npar  <- length(par_names)
ndisp <- length(disp_names)

result_table <- data.frame(
  Parameter  = c(par_names, disp_names),
  Estimate   = c(unlist(fixed_est), unlist(disp_est)),
  SE         = c(md0$fixed_sd, rep(NA, ndisp)),
  stringsAsFactors = FALSE
) %>%
  mutate(
    `$z$-value` = Estimate / SE,
    `$p$-value` = 2 * pnorm(abs(`$z$-value`), lower.tail = FALSE)
  )
result_table

latex_table <- xtable(result_table, type = "latex", align = "cccccc")
digits(latex_table) <- 3
print(latex_table,
      file                    = here::here("Data Analysis/Results/Table.tex"),
      include.rownames        = FALSE,
      sanitize.text.function  = identity,
      hline.after             = c(-1, 0, npar + ndisp))

# Random-effect correlation matrix
round(solve(md0$invSIGMA), 3)

saveRDS(md0, "Data Analysis/Results/Results.RDS")

# =============================================================================
# 10. Ti TABLE (estimated change-point per patient)
# =============================================================================

Ti_table <- data_cd4 %>%
  left_join(md0$Bi, by = "PATIENT") %>%
  left_join(md0$Ti, by = "PATIENT") %>%
  mutate(fitted_CD4 = fixed_est$alpha0 + disp_est$Lambda.alpha0 * a0i +
           fixed_est$alpha1 * time_CD4 + fixed_est$alpha2 * time_CD4^2) %>%
  group_by(PATIENT) %>%
  mutate(max_CD4 = max(fitted_CD4)) %>%
  dplyr::select(-c(time_CD4, y_CD4, fitted_CD4, y_logCD4)) %>%
  distinct() %>%
  mutate(Ti_est = exp(fixed_est$omega1 * tau2i + fixed_est$omega2 * max_CD4)) %>% 
  as.data.frame()
Ti_table
