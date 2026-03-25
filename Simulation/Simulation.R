# =============================================================================
# HIV/AIDS Joint Modeling: Simulation Study
# =============================================================================

library(tidyverse)
library(nlme)
library(optimx)
library(mvtnorm)
library(xtable)
library(glue)

setwd(here::here())

# Source all helper functions from the R/ directory
file.sources <- list.files("R", pattern = "\\.R$", full.names = TRUE)
sapply(file.sources, source, .GlobalEnv)


# =============================================================================
# 1. SIMULATION SETTINGS  ← edit here to change any simulation configuration
# =============================================================================

cfg <- list(

  # --- Run settings ---
  num_sim     = 100,     # number of simulation replicates
  num_patient = 100,     # number of patients per replicate
  censor_value = log10(40),  # left-censoring limit (log10 scale)
  out_dir     = here::here("Simulation", "SimResults I"),

  # --- Observation time points ---
  time = c(0.03, 0.27, 0.47, 0.93, 1.25, 1.43, 1.87, 2.17, 2.25, 2.56,
           2.80, 3.00, 3.15, 3.34, 3.40, 3.68, 3.87, 4.00, 4.20, 4.37,
           4.65, 4.78, 4.94, 5.04, 5.20, 5.45, 6.00, 6.10, 6.31, 6.53,
           6.87, 7.10, 7.28, 7.45, 7.76, 8.10, 8.42, 8.80, 9.10),

  # --- True fixed effects ---
  # Viral decay
  eta1 = 0.8, eta2 = 0.7, eta3 = 3,
  # Viral rebound
  beta1 = 1.8, beta2 = 1.0, beta3 = 3.0, beta4 = 1.2,
  # CD4
  alpha0 = 3.0, alpha1 = 0.2, alpha2 = -0.03,
  # Transition
  omega1 = -0.2, omega2 = 0.1,

  # --- True dispersion parameters ---
  sigma1 = 0.06,   # viral decay residual SD
  sigma2 = 0.10,   # viral rebound residual SD
  sigma3 = 0.40,   # CD4 residual SD
  sigma4 = 0.10,   # transition residual SD

  # --- Random-effect SDs (0 = not included) ---
  Lambda.eta   = c(0.1, 0.05, 0.0),   # eta1, eta2, eta3
  Lambda.beta  = c(0.1, 0.0, 0.1, 0.0),  # beta1, beta2, beta3, beta4
  Lambda.alpha = c(0.1, 0.0, 0.0),    # alpha0, alpha1, alpha2

  # --- Random-effect correlation (between non-zero effects) ---
  # Off-diagonal entry for Sigma[2, 3] (tau2i <-> b1i)
  Sigma_cor_23 = 0.8
)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 2. DERIVED OBJECTS FROM CONFIG
# =============================================================================

# Named Lambda vectors
names(cfg$Lambda.eta)   <- paste0("Lambda.eta",   seq_along(cfg$Lambda.eta))
names(cfg$Lambda.beta)  <- paste0("Lambda.beta",  seq_along(cfg$Lambda.beta))
names(cfg$Lambda.alpha) <- paste0("Lambda.alpha", seq(0, length(cfg$Lambda.alpha) - 1))

Lambda.all     <- c(cfg$Lambda.eta, cfg$Lambda.beta, cfg$Lambda.alpha)
Lambda.nonzero <- Lambda.all[Lambda.all != 0]
Lambda.zero    <- Lambda.all[Lambda.all == 0]
pos.nonzero    <- which(Lambda.all != 0)
n_randef       <- length(Lambda.nonzero)

# Random-effect covariance matrix (spherical base + one off-diagonal correlation)
Sigma <- diag(n_randef)
Sigma[2, 3] <- Sigma[3, 2] <- cfg$Sigma_cor_23

# Embed Sigma into full (10 x 10) matrix aligned with all random effects
newSigma <- diag(10)
newSigma[pos.nonzero, pos.nonzero] <- Sigma

# Named true-value lists (used as starting values and for results table)
par_names  <- c("eta1", "eta2", "eta3",
                "beta1", "beta2", "beta3", "beta4",
                "alpha0", "alpha1", "alpha2",
                "omega1", "omega2")
true_value <- c(cfg$eta1, cfg$eta2, cfg$eta3,
                cfg$beta1, cfg$beta2, cfg$beta3, cfg$beta4,
                cfg$alpha0, cfg$alpha1, cfg$alpha2,
                cfg$omega1, cfg$omega2)
true_value_par <- setNames(as.list(true_value), par_names)

# Dispersion parameter names (off-diagonal Sigma entries + Lambdas + sigmas)
upper_Sigma      <- Sigma
upper_Sigma[lower.tri(upper_Sigma, diag = TRUE)] <- 0
nonzero_offdiag  <- which(upper_Sigma != 0, arr.ind = TRUE)
Sigma_names      <- apply(nonzero_offdiag, 1,
                           function(idx) paste0("Sigma", idx[1], idx[2]))

disppar_names <- c(Sigma_names, names(Lambda.nonzero),
                   "y_decay_sigma", "y_rebound_sigma", "CD4_sigma", "time_sigma")
disp_value    <- c(Sigma[nonzero_offdiag], Lambda.nonzero,
                   cfg$sigma1, cfg$sigma2, cfg$sigma3, cfg$sigma4)

disp_value_par <- c(as.list(disp_value), as.list(Lambda.zero))
names(disp_value_par)[seq_along(disppar_names)] <- disppar_names


# =============================================================================
# 3. MODEL SPECIFICATION OBJECTS
# =============================================================================

Object_decay <- list(
  modeltype      = "nlme",
  response       = "y_decay",
  reg_equation   = "(eta1 + Lambda.eta1 * tau1i) + eta3 * exp(-(eta2 + Lambda.eta2 * tau2i) * time_decay)",
  distribution   = "normal",
  random_effects = c("tau1i", "tau2i"),
  fixed_param    = list(names        = c("eta1", "eta2", "eta3"),
                        start_values = c(true_value_par$eta1, true_value_par$eta2, true_value_par$eta3),
                        lower        = rep(-Inf, 3), upper = rep(Inf, 3)),
  disp_param     = list(names        = c("Lambda.eta1", "Lambda.eta2"),
                        start_values = c(cfg$Lambda.eta[1], cfg$Lambda.eta[2]),
                        lower        = c(0, 0), upper = c(Inf, Inf)),
  sigma          = disp_value_par$y_decay_sigma,
  left_censoring = list(indicator   = "censor_decay",
                        limit_value = cfg$censor_value,
                        method      = "tobit")
)

Object_rebound <- list(
  modeltype      = "nlme",
  response       = "y_rebound",
  reg_equation   = "(beta1 + Lambda.beta1 * b1i) * time_rebound / (time_rebound + exp(beta2 - (beta3 + Lambda.beta3 * b3i) * time_rebound)) + beta4",
  distribution   = "normal",
  random_effects = c("b1i", "b3i"),
  fixed_param    = list(names        = c("beta1", "beta2", "beta3", "beta4"),
                        start_values = c(true_value_par$beta1, true_value_par$beta2,
                                         true_value_par$beta3, true_value_par$beta4),
                        lower        = rep(-Inf, 4), upper = rep(Inf, 4)),
  disp_param     = list(names        = c("Lambda.beta1", "Lambda.beta3"),
                        start_values = c(cfg$Lambda.beta[1], cfg$Lambda.beta[3]),
                        lower        = c(0, 0), upper = c(Inf, Inf)),
  sigma          = disp_value_par$y_rebound_sigma,
  left_censoring = list(indicator   = "censor_rebound",
                        limit_value = cfg$censor_value,
                        method      = "tobit")
)

Object_CD4 <- list(
  modeltype      = "nlme",
  response       = "CD4",
  reg_equation   = "(alpha0 + Lambda.alpha0 * a0i) + alpha1 * time_CD4 + alpha2 * time_CD4^2",
  distribution   = "normal",
  random_effects = "a0i",
  fixed_param    = list(names        = paste0("alpha", 0:2),
                        start_values = c(true_value_par$alpha0, true_value_par$alpha1,
                                         true_value_par$alpha2),
                        lower        = rep(-Inf, 3), upper = rep(Inf, 3)),
  disp_param     = list(names        = "Lambda.alpha0",
                        start_values = cfg$Lambda.alpha[1],
                        lower        = 0, upper        = Inf),
  sigma          = disp_value_par$CD4_sigma
)

Object_trans <- list(
  modeltype    = "survival",
  response     = "time",
  event        = "trans",
  reg_equation = "omega1 * tau2i + omega2 * maxCD4",
  distribution = "log-normal",
  fixed_param  = list(names        = c("omega1", "omega2"),
                      start_values = c(true_value_par$omega1, true_value_par$omega2),
                      lower        = rep(-Inf, 2), upper = rep(Inf, 2)),
  sigma        = disp_value_par$time_sigma,
  truncated    = list(lower = "treatment_stop", upper = "first_reboundobs")
)

memList <- list(Object_decay, Object_rebound, Object_CD4, Object_trans)


# =============================================================================
# 4. SIMULATION LOOP
# =============================================================================

par_length     <- length(true_value)
disppar_length <- length(disp_value)

table_par   <- data.frame(matrix(NA, nrow = cfg$num_sim, ncol = par_length))
table_parsd <- data.frame(matrix(NA, nrow = cfg$num_sim, ncol = par_length))
table_disp  <- data.frame(matrix(NA, nrow = cfg$num_sim, ncol = disppar_length))
colnames(table_par) <- colnames(table_parsd) <- par_names
colnames(table_disp) <- disppar_names

i  <- 1
NC <- 0   # non-convergence counter

while (i <= cfg$num_sim) {

  cat(glue("Replicate {i} / {cfg$num_sim}\n"))

  # Simulate data
  data_sim <- create_simdata(cfg$time, cfg$num_patient, cfg$censor_value,
                             true_value_par, disp_value_par, newSigma, trt = TRUE)

  dataList <- list(data_sim$decay, data_sim$rebound, data_sim$CD4, data_sim$trans)

  # Fit joint model
  md0 <- try(
    fit_multi_mems(memList      = memList,
                   dataList     = dataList,
                   subject_id   = "ID",
                   randeff_info = list(distribution = "normal", degree = 0),
                   cov_method   = "spherical",
                   loglike_tol  = 1e-3, par_tol = 1e-3,
                   Silent       = TRUE, iterMax = 30,
                   REML         = FALSE, naive = FALSE, adjust = TRUE),
    silent = TRUE
  )

  # Skip replicate on convergence failure or missing SEs
  if (inherits(md0, "try-error") || any(is.na(md0$fixed_sd))) {
    NC <- NC + 1
    cat(glue("  Non-convergence (total so far: {NC})\n"))
    next
  }

  # Store estimates
  table_par[i, ]   <- round(md0$fixed_est, 4)
  table_parsd[i, ] <- round(as.numeric(md0$fixed_sd), 4)
  table_disp[i, ]  <- c(solve(md0$invSIGMA)[2, 3], round(md0$disp_est, 4))

  # Checkpoint save after every accepted replicate
  saveRDS(table_par,   file.path(cfg$out_dir, "table_par.RDS"))
  saveRDS(table_parsd, file.path(cfg$out_dir, "table_parsd.RDS"))
  saveRDS(table_disp,  file.path(cfg$out_dir, "table_disp.RDS"))
  saveRDS(md0,         file.path(cfg$out_dir, "example_result.RDS"))

  i <- i + 1
}

cat(glue("\nSimulation complete. Non-convergences: {NC} / {cfg$num_sim + NC}\n"))


# =============================================================================
# 5. RESULTS TABLE
# =============================================================================

table_par   <- readRDS(file.path(cfg$out_dir, "table_par.RDS"))
table_parsd <- readRDS(file.path(cfg$out_dir, "table_parsd.RDS"))
table_disp  <- readRDS(file.path(cfg$out_dir, "table_disp.RDS"))

simresult <- simtable(table_par, table_parsd, par_names, true_value,
                      table_disp, disppar_names, disp_value)

simresult$Parameter <- c(
  "$\\eta_1$", "$\\eta_2$", "$\\eta_3$",
  "$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$",
  "$\\alpha_0$", "$\\alpha_1$", "$\\alpha_2$",
  "$\\omega_1$", "$\\omega_2$",
  "$\\Sigma_{2,4}$",
  "$\\Lambda_{\\tau11}$", "$\\Lambda_{\\tau22}$",
  "$\\Lambda_{b11}$", "$\\Lambda_{b33}$",
  "$\\Lambda_{a11}$",
  "$\\sigma_1$", "$\\sigma_2$", "$\\sigma_3$", "$\\sigma_4$"
)

latex_table <- xtable(simresult, type = "latex", align = "ccccccccc")
digits(latex_table) <- c(0, 2, 2, 3, 3, 3, 3, 3, 2)
print(latex_table,
      file                   = file.path(cfg$out_dir, "latex_table.tex"),
      include.rownames       = FALSE,
      sanitize.text.function = identity,
      hline.after            = c(-1, 0, length(true_value), nrow(simresult)))
