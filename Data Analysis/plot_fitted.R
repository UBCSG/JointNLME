# =============================================================================
# HIV/AIDS Joint Modeling: Individual Fitted Trajectory Plots
# =============================================================================

library(tidyverse)
library(berryFunctions)  

outdir <- here::here("Data Analysis/Results/FittedPlots")
dir.create(outdir, showWarnings = FALSE)

data_decay <- read.csv(here::here("Data", "cleaned_data_decay.csv"))
md0 <- readRDS("Data Analysis/Results/Results.RDS")
fixed_est <- as.list(md0$fixed_est)
disp_est  <- as.list(md0$disp_est)

# =============================================================================
# HELPER: fitted viral-load curves
# =============================================================================

# Decay-phase curve (ART on, t in [0, Ti_est])
curve_decay <- function(t, fe, de, re) {
  (fe$eta1 + de$Lambda.eta1 * re$tau1i) +
    fe$eta3 * exp(-(fe$eta2 + de$Lambda.eta2 * re$tau2i) * t)
}

# Rebound-phase curve (ART off, t shifted so origin = Ti_est)
curve_rebound <- function(t, Ti, fe, de, re) {
  s <- t - Ti   # time since change-point
  (fe$beta1 + de$Lambda.beta1 * re$b1i) * s /
    (s + exp(fe$beta2 - (fe$beta3 + de$Lambda.beta3 * re$b3i) * s)) +
    fe$beta4 / (1 + exp(fe$beta5 * s))
}


# =============================================================================
# HELPER: standard individual plot (all patients)
# =============================================================================

plot_annotated <- function(patient_id, data_decay, md0, fixed_est, disp_est,
                           x_offset  = 0.03,   # point nudge (years) after Ti_est
                           zoom_half = 0.25,   # half-width of zoom window (years)
                           base_size = 20) {
  
  fe <- fixed_est
  de <- disp_est
  re <- dplyr::filter(md0$Bi, PATIENT == patient_id)
  
  subdat         <- data_decay %>% filter(PATIENT == patient_id)
  treatment_stop <- unique(subdat$treatment_stop)
  Ti_est         <- md0$Ti %>% filter(PATIENT == patient_id) %>%
    pull(Ti0_sim) + treatment_stop
  upper_bound    <- unique(subdat$first_reboundobs)
  top            <- 1.1 * max(subdat$y_decay_halfDL)
  ymax           <- max(subdat$y_decay_halfDL)
  
  # Offset rebound points so they do not overlap the Ti_est vline
  subdat <- subdat %>%
    mutate(time_plot = ifelse(time_decay > Ti_est, time_decay + x_offset, time_decay))
  t_end <- max(subdat$time_plot)
  
  # Flip annotation side when a landmark sits near a plot edge (avoids clipping)
  Ti_label_x <- ifelse(Ti_est / t_end > 0.85, 0.96 * Ti_est, 1.04 * Ti_est)
  Ti_hjust   <- ifelse(Ti_est / t_end > 0.85, 1.1, -0.1)
  ub_label_x <- ifelse(upper_bound / t_end > 0.85, 0.92 * upper_bound, 1.08 * upper_bound)
  ts_label_x <- ifelse(treatment_stop / t_end < 0.15, 1.07 * treatment_stop, 0.93 * treatment_stop)
  
  p <- ggplot() +
    theme_classic(base_size = base_size) +
    xlim(0, t_end) + ylim(0, top) +
    
    # Estimated change-point (Ti_est)
    geom_vline(xintercept = Ti_est, linetype = 2) +
    annotate("text",
             x = Ti_label_x, y = 0.85 * top, hjust = Ti_hjust, size = 6,
             label = expression(hat(T)[JM]), parse = TRUE) +
    geom_segment(aes(x = Ti_label_x, y = 0.80 * top,
                     xend = Ti_est,  yend = 0.75 * top),
                 arrow = arrow(length = unit(0.3, "cm"))) +
    
    # Observed data
    geom_point(data = subdat,
               aes(x = time_plot, y = y_decay_halfDL, shape = factor(censor_decay)),
               size = 2) +
    scale_shape_manual(name   = "Censored",
                       values = c("0" = 19, "1" = 1),
                       labels = c("0" = "No", "1" = "Yes")) +
    
    # Fitted decay and rebound curves
    geom_function(fun  = \(x) curve_decay(x, fe, de, re),
                  xlim = c(0, Ti_est)) +
    geom_function(fun  = \(x) curve_rebound(x - x_offset, Ti_est, fe, de, re),
                  xlim = c(Ti_est + x_offset, t_end)) +
    
    # Smooth visual bridge across the x_offset gap
    geom_curve(aes(x = Ti_est, y = 0.5,
                   xend = Ti_est + x_offset, yend = log10(40) + 0.05),
               curvature = 0.08) +
    
    # Detection limit
    geom_hline(yintercept = log10(40), linetype = "dotted") +
    annotate("text", x = 0.9 * t_end, y = 1.4, label = "Detection limit", size = 6) +
    
    # Latest change-point (upper_bound)
    geom_vline(xintercept = upper_bound, linetype = "dotted") +
    annotate("text",
             x = ub_label_x, y = 0.63 * ymax,
             label = "Latest\nchange\npoint", size = 6) +
    geom_segment(aes(x    = ub_label_x,         y    = 0.48 * ymax,
                     xend = 1.01 * upper_bound,  yend = 0.40 * ymax),
                 arrow = arrow()) +
    
    # Treatment stop
    geom_vline(xintercept = treatment_stop, linetype = "dotted") +
    annotate("text",
             x = ts_label_x, y = 0.63 * ymax,
             label = "Treatment\nstop\ntime", size = 6) +
    geom_segment(aes(x    = ts_label_x,            y    = 0.48 * ymax,
                     xend = 0.99 * treatment_stop,  yend = 0.40 * ymax),
                 arrow = arrow()) +
    
    xlab("Time in years") +
    ylab(bquote("Viral load (in" ~ log[10] ~ "-scale)"))
  
  p
}


# =============================================================================
# GENERATE PLOTS — annotated dual-panel for every patient
# =============================================================================

for (i in unique(data_decay$PATIENT)) {
  p <- plot_annotated(i, data_decay, md0, fixed_est, disp_est)
  ggsave(file.path(outdir, sprintf("ID%d.png", i)), p,
         width = 20, height = 6)
  cat("Saved patient", i, "\n")
}
