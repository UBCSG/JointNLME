# ==========================================================================
# HIV/AIDS Joint Modeling: Visualization
# ==========================================================================

library(tidyverse)
library(patchwork)
library(scales)
library(here)

# 1. Data Loading ---------------------------------------------------------

# Load cleaned datasets from the Data directory
data <- read.csv(here::here("Data", "cleaned_data.csv"))
data_cd4 <- read.csv(here::here("Data", "cleaned_data_cd4.csv"))
data_decay <- read.csv(here::here("Data", "cleaned_data_decay.csv"))

# 2. Individual Patient Plots ---------------------------------------------

# Get list of unique patients to iterate through
uniqueID <- unique(data$PATIENT)

for (i in uniqueID) {
  # Filter data for current subject
  subdat <- data %>% filter(PATIENT == i)
  subdat_cd4 <- data_cd4 %>% filter(PATIENT == i)
  
  # Rescale CD4 values to the Viral Load range for secondary axis visualization
  # This allows both biomarkers to be plotted on the same physical space
  subdat_cd4$y_cd4_recale <- scales::rescale(
    subdat_cd4$y_logCD4, 
    to = c(min(subdat$log_rna_halfDL), max(subdat$log_rna_halfDL))
  )
  
  # Define plot boundary constants
  t_end <- max(subdat$time_years)
  Ti <- unique(subdat$treatment_stop_years)
  upper_bound <- unique(subdat$first_reboundobs)
  top <- 1.1 * max(subdat$log_rna_halfDL)
  treatment_stop <- unique(subdat$treatment_stop_years)
  
  # Identify peak CD4 timing between treatment stop and viral rebound
  subdat_cd4_window <- subdat_cd4 %>% 
    filter(time_CD4 <= upper_bound & time_CD4 >= treatment_stop)
  
  t_max_cd4 <- max(subdat_cd4_window[which(subdat_cd4_window$y_logCD4 == max(subdat_cd4_window$y_logCD4)), ]$time_CD4)
  
  # Add labels for legend grouping
  subdat$measure <- "Viral load"
  subdat_cd4$measure <- "CD4"
  
  # Construct individual plot
  p <- ggplot() +
    xlim(0, t_end) +
    ylim(0, top) +
    theme_classic() +
    theme(text = element_text(size = 16)) +
    
    # Layer 1: Data points (Viral Load circle/open circle, CD4 triangle)
    geom_point(data = subdat, aes(x = time_years, y = log_rna_halfDL, shape = factor(censor)), size = 3) +
    geom_point(data = subdat_cd4, aes(x = time_CD4, y = y_cd4_recale, shape = measure), size = 3) +
    
    # Layer 2: Connecting lines
    geom_line(data = subdat, aes(x = time_years, y = log_rna_halfDL)) +
    geom_line(data = subdat_cd4, aes(x = time_CD4, y = y_cd4_recale)) +
    
    # Layer 3: Annotations and vertical markers
    geom_hline(yintercept = log10(40), linetype = 'dashed') + # Detection limit
    geom_text(data = subdat, x = 0.9 * t_end, y = 1.4, size = 6, label = "detection limit") +
    
    geom_vline(xintercept = upper_bound, linetype = 'dotted') +
    geom_text(data = subdat, x = 1.05 * upper_bound, y = 2.5, size = 6, label = "Latest\nchange\npoint") +
    
    geom_vline(xintercept = Ti, linetype = 'dotted') +
    geom_text(data = subdat, x = 0.95 * Ti, y = 2.5, size = 6, label = "Treatment\nstop\ntime") +
    
    # Manual aesthetic controls
    scale_shape_manual(
      name = "Measurement",
      values = c("0" = 19, "1" = 1, "CD4" = 17),
      labels = c("0" = "Viral load", "1" = "Viral load (censored)", "CD4" = "CD4")
    ) +
    xlab("Time in years") +
    scale_y_continuous(
      name = bquote("Viral load (in" ~ log[10] ~ "-scale)"),
      # Secondary axis handles the back-transformation for display labels
      sec.axis = sec_axis(
        ~ scales::rescale(., 
                          from = c(min(subdat$log_rna_halfDL), max(subdat$log_rna_halfDL)), 
                          to = c(min(subdat_cd4$y_logCD4), max(subdat_cd4$y_logCD4))),
        name = bquote("CD4 count (in" ~ log[10] ~ "-scale)")
      )
    )
  
  # Save individual plot
  ggsave(paste0("Figures/Ind_plot/ID", i, ".png"), width = 20, height = 6)
  cat("Processed Patient ID:", i, "\n")
}

# 3. Population Viral Load Trajectories -----------------------------------

# Spaghetti plot of viral load for all patients
ggplot(data_decay, aes(x = time_decay, y = y_decay_halfDL)) + 
  geom_point(aes(shape = factor(censor_decay)), size = 2) +
  geom_line(aes(group = PATIENT)) +
  geom_hline(yintercept = log10(40), linetype = "dotted") + 
  annotate("text", x = 5, y = 1.8, label = "Detection Limit", size = 5) +
  scale_x_continuous("Time (in years)") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10] ~ "-scale)")) +
  scale_shape_manual(name = "Censored", values = c(16, 1), labels = c('No', 'Yes')) +
  theme_classic() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"))

ggsave("Figures/Viral load all patients.png", width = 10, height = 5)

# 4. Highlighted Trajectories (3 Subjects) --------------------------------

# Create annotation data for specific treatment stop labels
ann_data <- data.frame(
  PATIENT = c(1, 8, 11),
  label = c("Stop Time 1", "Stop Time 2", "Stop Time 3"),
  treatment_stop = c(1.643836, 3.767123, 1.586301),
  x_pos = c(2, 3.4, 1.2), 
  y_pos = c(7, 7, 7)         
)

# Identify endpoint positions for patient bold labels
ann_data2 <- data_decay %>%
  filter(PATIENT %in% c(1, 8, 11)) %>%
  group_by(PATIENT) %>%
  filter(time_decay == max(time_decay)) %>%
  mutate(
    label = case_when(
      PATIENT == 1  ~ "Patient 1",
      PATIENT == 8  ~ "Patient 2",
      PATIENT == 11 ~ "Patient 3"
    ),
    x_pos = time_decay + 0.1,
    y_pos = y_decay_halfDL
  )

# Plot trajectories for the 3 selected patients with custom markers
ggplot(data_decay %>% filter(PATIENT %in% c(1, 11, 8)), aes(x = time_decay, y = y_decay_halfDL)) + 
  geom_point(aes(shape = factor(censor_decay)), size = 2) +
  geom_line(aes(group = PATIENT)) +
  # Add vertical lines for ART cessation
  geom_vline(aes(xintercept = treatment_stop), linetype = "dashed", color = "darkgrey") +
  # Add stop time labels and arrows
  geom_text(data = ann_data, aes(x = x_pos, y = y_pos, label = label), size = 4) +
  geom_curve(data = ann_data, 
             aes(x = x_pos + 0.1, y = y_pos - 0.2, xend = treatment_stop, yend = y_pos - 0.5),
             arrow = arrow(length = unit(0.2, "cm")), 
             curvature = -0.2) +
  # Label endpoints of trajectories
  geom_text(data = ann_data2, aes(x = x_pos - 0.5, y = y_pos + 0.3, label = label), 
            hjust = 0, size = 4, fontface = "bold") +
  
  geom_hline(yintercept = log10(40), linetype = "dotted") + 
  annotate("text", x = 4.5, y = 1.8, label = "Detection Limit", size = 5) +
  scale_x_continuous("Time (in years)") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10] ~ "-scale)")) +
  scale_shape_manual(name = "Censored", values = c(16, 1), labels = c('No', 'Yes')) +
  theme_classic() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_line(colour = "grey90"), 
        panel.grid.minor = element_line(colour = "grey90"))

ggsave("Figures/Viral load 3 patients.png", width = 10, height = 5)
