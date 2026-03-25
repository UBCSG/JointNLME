# =============================================================================
# HIV/AIDS Joint Modeling: Example Fitted Plot (Single Patient from Simulation)
# =============================================================================

library(tidyverse)
library(ggplot2)

# =============================================================================
# 1. LOAD RESULTS & EXTRACT PATIENT DATA
# =============================================================================

dir <- here::here("Simulation", "SimResults I")
md0 <- readRDS(here::here(dir, "example_result.RDS"))

q <- 1  # patient index to plot

fixed_est <- as.list(md0$fixed_est)
disp_est  <- as.list(md0$disp_est)
randef    <- as.list(md0$Bi %>% filter(ID == q))
censor_value <- log10(40)

subdata_decay   <- md0$dataList[[1]] %>% filter(ID == q)
subdata_rebound <- md0$dataList[[2]] %>% filter(ID == q)
subTidata       <- md0$dataList[[4]] %>% filter(ID == q)

trt_stop <- unique(subdata_decay$treatment_stop)
tq2      <- unique(subdata_rebound$first_reboundobs)  # latest change-point upper bound
t_end    <- max(subdata_decay$time_decay)

# Estimated change-point and 95% CI (on log-normal scale)
Ti0_est       <- md0$Ti$Ti0_sim[q]
Ti_est        <- trt_stop + Ti0_est
Ti_est_lower  <- max(trt_stop, trt_stop + Ti0_est * exp(-1.96 * disp_est$time_sigma))
Ti_est_upper  <- min(tq2,      trt_stop + Ti0_est * exp( 1.96 * disp_est$time_sigma))

# Half detection-limit substitution for censored observations
subdata_decay <- subdata_decay %>%
  mutate(y_plot = ifelse(censor_decay == 1, log10(0.5 * 10^y_decay), y_decay))

# x-value where decay curve ends and the bridge curve begins
x_decay_end <- 2.2


# =============================================================================
# 2. CURVE FUNCTIONS (patient-specific fitted values)
# =============================================================================

decay_curve <- function(x) {
  (fixed_est$eta1 + disp_est$Lambda.eta1 * randef$tau1i) +
    fixed_est$eta3 * exp(-(fixed_est$eta2 + disp_est$Lambda.eta2 * randef$tau2i) * x)
}

rebound_curve <- function(x) {
  s <- x - Ti_est
  (fixed_est$beta1 + disp_est$Lambda.beta1 * randef$b1i) * s /
    (s + exp(fixed_est$beta2 - (fixed_est$beta3 + disp_est$Lambda.beta3 * randef$b3i) * s)) +
    fixed_est$beta4
}


# =============================================================================
# 3. PLOT
# =============================================================================

ggplot(subdata_decay, aes(x = time_decay, y = y_plot)) +

  # Shaded CI band for Ti estimate
  annotate("rect",
           xmin = Ti_est_lower, xmax = Ti_est_upper, ymin = 0, ymax = 5,
           alpha = 0.3) +

  # Observed data points
  geom_point(aes(shape = factor(censor_decay)), size = 2) +
  scale_shape_manual(name   = "Censored",
                     values = c("0" = 19, "1" = 1),
                     labels = c("0" = "No", "1" = "Yes")) +

  # Fitted curves
  geom_function(fun = decay_curve,   xlim = c(0, x_decay_end)) +
  geom_function(fun = rebound_curve, xlim = c(Ti_est, t_end)) +

  # Smooth bridge connecting decay endpoint to rebound start
  geom_curve(aes(x    = x_decay_end,
                 y    = decay_curve(x_decay_end),
                 xend = Ti_est,
                 yend = fixed_est$beta4),
             curvature = 0.05) +

  # Detection limit
  geom_hline(yintercept = censor_value, linetype = "dotted") +
  annotate("text", x = 8, y = 1.1 * censor_value,
           label = "Detection~Limit", parse = TRUE, hjust = 1.1, size = 6) +

  # Vertical reference lines with annotations
  # — Latest change-point upper bound (tq2)
  geom_vline(xintercept = tq2, linetype = "dotted") +
  annotate("text", x = 1.06 * tq2, y = 3,
           label = "t[iq2]", parse = TRUE, hjust = 1.1, size = 6) +
  geom_segment(aes(x = 1.03 * tq2, y = 2.9, xend = 1.01 * tq2, yend = 2.7),
               arrow = arrow(length = unit(0.3, "cm"))) +

  # — Treatment stop (tq1)
  geom_vline(xintercept = trt_stop, linetype = "dotted") +
  annotate("text", x = 0.98 * trt_stop, y = 3,
           label = "t[iq1]", parse = TRUE, hjust = 1.1, size = 6) +
  geom_segment(aes(x = 0.97 * trt_stop, y = 2.8, xend = 0.99 * trt_stop, yend = 2.6),
               arrow = arrow(length = unit(0.3, "cm"))) +

  # — Estimated change-point (T_hat_JM)
  geom_vline(xintercept = Ti_est, linetype = 2) +
  annotate("text", x = 0.985 * Ti_est, y = 3,
           label = expression(hat(T)[JM]), parse = TRUE, hjust = 1.1, size = 6) +
  geom_segment(aes(x = 0.98 * Ti_est, y = 2.8, xend = 0.99 * Ti_est, yend = 2.6),
               arrow = arrow(length = unit(0.3, "cm"))) +

  # — True change-point (T_i)
  geom_vline(xintercept = subTidata$Ti_true, linetype = 2) +
  annotate("text", x = 1.05 * subTidata$Ti_true, y = 3,
           label = expression(T[i]), parse = TRUE, hjust = 1.1, size = 6) +
  geom_segment(aes(x = 1.03 * subTidata$Ti_true, y = 2.9, xend = 1.01 * subTidata$Ti_true, yend = 2.7),
               arrow = arrow(length = unit(0.3, "cm"))) +

  theme_classic(base_size = 20) +
  ylim(0, 5) +
  scale_x_continuous("Time") +
  ylab(bquote("Viral load (in" ~ log[10] ~ "-scale)"))

ggsave(here::here(dir, "example_fitted.png"), width = 20, height = 6)
