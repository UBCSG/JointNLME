library(tidyverse)
library(ggplot2)
library(ggpubr)
library (berryFunctions)
setwd(here::here())

setwd("R")
file.sources = list.files(pattern = "*.R")
sapply(file.sources, source, .GlobalEnv)


# Import the raw data
raw_data <- read.csv(here::here("Data","raw_data.csv"), header = TRUE)

# Remove NA's
data <- subset(raw_data, !is.na(days_from_seroco) & !is.na(RNA_V))

data <- data %>% 
  mutate(log_cd4 = ifelse(is.na(CD4_V), NA, log(CD4_V))) %>% 
  mutate(censor = ifelse(RNA_V <= 40, 1, 0)) %>% 
  mutate(log_rna = ifelse(censor == 1, log10(40), log10(RNA_V))) %>% 
  mutate(log_rna_halfDL = ifelse(RNA_V <= 40, log10(20), log10(RNA_V))) %>% 
  filter(days_from_seroco >= art_start) %>% 
  mutate(time_days = days_from_seroco - art_start) %>% 
  mutate(treatment_stop = ati_start - art_start) %>% 
  mutate(limit = 0) %>% 
  mutate(time_years = time_days/365) %>% 
  mutate(treatment_stop_years = treatment_stop / 365) %>% 
  arrange(PATIENT, time_days)

data_cd4 <- data %>% filter(!is.na(log_cd4))
max_cd4 <- data_cd4 %>% 
  group_by(PATIENT) %>% 
  filter(log_cd4 == max(log_cd4)) %>% 
  dplyr::select(PATIENT, days_from_seroco, log_cd4)
colnames(max_cd4)[2] <- "time_log_cd4_max"
colnames(max_cd4)[3] <- "log_cd4_max"

data <- merge(data, max_cd4) %>% 
  mutate(last_decayobs = 0)

# Manually remove some outliers
data <- data %>% 
  filter(!(PATIENT == 2 & time_days > 1400)) %>% 
  filter(!(PATIENT == 3 & time_days == 875)) %>% 
  filter(!(PATIENT == 3 & time_days > 1100)) %>% 
  filter(!(PATIENT == 4 & time_days == 793)) %>% 
  filter(!(PATIENT == 5 & time_days == 1040)) %>% 
  filter(!(PATIENT == 5 & time_days > 1100)) %>% 
  filter(!(PATIENT == 7 & time_days == 389)) %>% 
  filter(!(PATIENT == 7 & time_days > 900)) %>% 
  filter(!(PATIENT == 8 & time_days > 1800)) %>% 
  filter(!(PATIENT == 9 & time_days > 900)) %>% 
  filter(!(PATIENT == 11 & time_days == 670)) %>% 
  filter(!(PATIENT == 11 & time_days == 706)) %>% 
  filter(!(PATIENT == 11 & time_days > 1200)) %>% 
  filter(!(PATIENT == 12 & time_days == 732)) %>% 
  filter(!(PATIENT == 12 & time_days == 942)) %>% 
  filter(!(PATIENT == 12 & time_days > 1200)) %>% 
  filter(!(PATIENT == 14 & time_days == 784)) %>% 
  filter(!(PATIENT == 14 & time_days > 900)) %>% 
  filter(!(PATIENT == 15 & time_days > 1100)) %>% 
  filter(!(PATIENT == 16 & time_days == 630)) %>% 
  filter(!(PATIENT == 16 & time_days == 679)) %>% 
  filter(!(PATIENT == 16 & time_days > 1200)) %>% 
  filter(!(PATIENT == 17 & time_days == 698)) %>% 
  filter(!(PATIENT == 17 & time_days > 1000)) %>% 
  filter(!(PATIENT == 18 & time_days > 1500)) %>% 
  filter(!(PATIENT == 19 & time_days > 800)) %>% 
  filter(!(PATIENT == 20 & time_days > 1200)) %>% 
  filter(!(PATIENT == 21 & time_days == 550)) %>% 
  filter(!(PATIENT == 21 & time_days == 579)) %>% 
  filter(!(PATIENT == 21 & time_days == 818)) %>% 
  filter(!(PATIENT == 21 & time_days == 865)) %>% 
  filter(!(PATIENT == 21 & time_days > 750)) %>% 
  filter(!(PATIENT == 23 & time_days > 450)) %>% 
  filter(!(PATIENT == 24 & time_days > 1000)) %>% 
  filter(!(PATIENT == 24 & time_days == 454)) %>% 
  filter(!(PATIENT == 25 & time_days > 1000)) %>% 
  filter(!(PATIENT == 26)) %>% 
  filter(!(PATIENT == 27)) %>% 
  filter(!(PATIENT == 28 & time_days > 1100)) %>% 
  filter(!(PATIENT == 29 & time_days == 658)) %>% 
  filter(!(PATIENT == 29 & time_days > 1000)) %>% 
  filter(!(PATIENT == 30 & time_days > 1200)) %>% 
  filter(!(PATIENT == 32 & time_days == 503)) %>% 
  filter(!(PATIENT == 32 & time_days == 559)) %>% 
  filter(!(PATIENT == 33 & time_days > 1000)) %>% 
  filter(!(PATIENT == 35 & time_days > 1000)) %>% 
  filter(!(PATIENT == 36 & time_days > 800)) %>% 
  filter(!(PATIENT == 37 & time_days > 1000)) %>% 
  filter(!(PATIENT == 38 & time_days > 1200)) %>% 
  filter(!(PATIENT == 40 & time_days > 1200)) %>% 
  filter(!(PATIENT == 41)) %>% 
  filter(!(PATIENT == 42 & time_days > 1000)) %>% 
  filter(!(PATIENT == 43 & time_days > 1200)) %>% 
  filter(!(PATIENT == 46 & time_days == 749)) %>% 
  filter(!(PATIENT == 47 & time_days > 800)) %>% 
  filter(!(PATIENT == 50 & time_days > 1500)) %>% 
  filter(!(PATIENT == 52 & time_days > 900)) %>% 
  filter(!(PATIENT == 53 & time_days == 445)) %>% 
  filter(!(PATIENT == 53 & time_days > 1000)) %>% 
  filter(!(PATIENT == 56 & time_days == 445)) %>% 
  filter(!(PATIENT == 56 & time_days > 1000)) %>% 
  filter(!(PATIENT == 57 & time_days > 1500)) %>% 
  filter(!(PATIENT == 57 & time_days == 982)) %>% 
  filter(!(PATIENT == 61 & time_days > 1100)) %>%
  filter(!(PATIENT == 62)) %>%
  filter(!(PATIENT == 64 & time_days > 1000)) %>%
  filter(!(PATIENT == 66 & time_days > 1200)) %>%
  filter(!(PATIENT == 67 & time_days > 800)) %>%
  filter(!(PATIENT == 69 & time_days > 800)) %>%
  filter(!(PATIENT == 70 & time_days > 800)) %>%
  filter(!(PATIENT == 71)) %>%
  filter(!(PATIENT == 73 & time_days == 210)) %>%
  filter(!(PATIENT == 73 & time_days == 232)) %>%
  filter(!(PATIENT == 73 & time_days > 1000)) %>%
  filter(!(PATIENT == 75 & time_days == 623)) %>%
  filter(!(PATIENT == 75 & time_days > 1200)) %>%
  filter(!(PATIENT == 76 & time_days > 1200)) %>%
  filter(!(PATIENT == 76 & time_days == 786)) 
  # mutate(treatment_stop = ifelse(PATIENT == 26, 292, treatment_stop)) %>% 
  # mutate(treatment_stop = ifelse(PATIENT == 41, 227, treatment_stop)) %>% 
  # mutate(treatment_stop = ifelse(PATIENT == 62, 243, treatment_stop)) %>% 
  # mutate(treatment_stop = ifelse(PATIENT == 71, 591, treatment_stop)) %>% 
  # mutate(treatment = ifelse(time_days <= treatment_stop, 1, 0))

# last_decayobs <- data %>%
#   filter(treatment == 1, !censor) %>%   # keep uncensored and above detection
#   group_by(PATIENT) %>%
#   slice_max(time_years) %>%
#   ungroup() %>% 
#   dplyr::select(PATIENT, time_years) %>% 
#   mutate(time_years = ifelse(time_years == 0, 0.1, time_years))
# colnames(last_decayobs)[2] <- "last_decayobs"
# data <- merge(data, last_decayobs) %>% 
#   mutate(last_decayobs = ifelse(PATIENT == 26, 104/365, last_decayobs)) %>% 
#   mutate(last_decayobs = ifelse(PATIENT == 33, 95/365, last_decayobs)) %>% 
#   mutate(last_decayobs = ifelse(PATIENT == 41, 93/365, last_decayobs)) %>% 
#   mutate(last_decayobs = ifelse(PATIENT == 49, 29/365, last_decayobs)) %>% 
#   mutate(last_decayobs = ifelse(PATIENT == 62, 159/365, last_decayobs)) %>% 
#   mutate(last_decayobs = ifelse(PATIENT == 71, 108/365, last_decayobs)) %>% 
#   mutate(last_decayobs = ifelse(PATIENT == 75, 122/365, last_decayobs)) 

first_reboundobs <- data %>%
  filter(treatment == 0, !censor) %>%   # keep uncensored and above detection
  group_by(PATIENT) %>%
  slice_min(time_years) %>%
  ungroup() %>% 
  dplyr::select(PATIENT, time_years)
colnames(first_reboundobs)[2] <- "first_reboundobs"
data <- merge(data, first_reboundobs) %>% 
  mutate(first_reboundobs = ifelse(PATIENT == 11, 782/365, first_reboundobs)) %>% 
  mutate(first_reboundobs = ifelse(PATIENT == 18, 752/365, first_reboundobs)) %>% 
  mutate(first_reboundobs = ifelse(PATIENT == 26, 348/365, first_reboundobs)) %>% 
  mutate(first_reboundobs = ifelse(PATIENT == 32, 1042/365, first_reboundobs)) %>% 
  mutate(first_reboundobs = ifelse(PATIENT == 41, 296/365, first_reboundobs)) %>% 
  mutate(first_reboundobs = ifelse(PATIENT == 56, 581/365, first_reboundobs)) %>% 
  mutate(first_reboundobs = ifelse(PATIENT == 62, 327/365, first_reboundobs)) %>% 
  mutate(first_reboundobs = ifelse(PATIENT == 71, 624/365, first_reboundobs)) %>% 
  arrange(PATIENT, time_years)

write.csv(data, here::here("Data", "cleaned_data.csv"), row.names = FALSE)
data <- read.csv(here::here("Data", "cleaned_data.csv"))
data_cd4 <- read.csv(here::here("Data", "cleaned_data_cd4.csv"))

uniqueID <- unique(data$PATIENT)
for(i in uniqueID){
  subdat <- data %>% filter(PATIENT == i)
  subdat_cd4 <- dplyr::filter(data_cd4, PATIENT == i)
  subdat_cd4$y_cd4_recale <- rescale(subdat_cd4$y_logCD4, 
                                     from = c(min(subdat_cd4$y_logCD4), max(subdat_cd4$y_logCD4)), 
                                     to = c(min(subdat$log_rna_halfDL), max(subdat$log_rna_halfDL)))
  t_end <- max(subdat$time_years)
  Ti <- unique(subdat$treatment_stop_years)
  upper_bound <- unique(subdat$first_reboundobs)
  subdat_2 <- subdat %>% filter(censor == 0)
  y <- median(subdat_2$log_rna)
  top <- 1.1*max(subdat$log_rna_halfDL)
  treatment_stop <- unique(subdat$treatment_stop_years)
  subdat_cd4_2 <- subdat_cd4 %>% filter(time_CD4 <= upper_bound & time_CD4 >= treatment_stop)
  t_max_cd4 <- max(subdat_cd4_2[which(subdat_cd4_2$y_logCD4 == max(subdat_cd4_2$y_logCD4)), ]$time_CD4)
  # p <- ggplot() +
  #   xlim(0, t_end) +
  #   ylim(1, 8) + ggtitle(paste("Paitient ID = ", i)) +
  #   # change background colors
  #   geom_rect(aes(ymin=1, ymax=8, xmin=0, xmax=Ti, fill="during ART"), alpha = .8) +
  #   geom_rect(aes(ymin=1, ymax=8, xmin=Ti, xmax = t_end, fill="after ART stop"), alpha = .8) +
  #   scale_fill_brewer(palette = 'Pastel2', name = 'Treatment Stage') +
  #   theme_classic() +
  #   # plot data
  #   geom_point(data = subdat, aes(x = time_years, y = log_rna_halfDL, colour = "log10(RNA)", shape = factor(censor))) +
  #   geom_point(data = subdat_cd4, aes(x = time_CD4, y = y_logCD4, colour = "log(CD4)")) +
  #   geom_line(data = subdat, aes(x = time_years, y = log_rna_halfDL, colour = "log10(RNA)")) +
  #   geom_line(data = subdat_cd4, aes(x = time_CD4, y = y_logCD4, colour = "log(CD4)")) +
  #   geom_hline(yintercept = log10(40), linetype = 'dashed') +
  #   geom_text(data = subdat, x = 9/10*t_end, y = 1.4, label = "detection limit") +
  #   geom_vline(xintercept = lower_bound, linetype = 'dashed') +
  #   geom_text(data = subdat, x = 1.1*lower_bound, y = 1.2*first_obs, label = "lower bound") +
  #   geom_vline(xintercept = upper_bound, linetype = 'dashed') +
  #   geom_text(data = subdat, x = 1.1*upper_bound, y = 1.2*first_obs, label = "upper bound") +
  #   scale_shape_manual(name = "Censor", values=c(19, 1), labels = c("No", "Yes")) +
  #   scale_colour_manual(name = "Biomarker",
  #                       values = c("blue", "black")) +
  #   xlab("Time in years") + ylab("Value")
  p <- ggplot() +
    xlim(0, t_end) +
    ylim(0, top) + ggtitle(paste("Paitient ID = ", i)) +
    # change background colors
    geom_rect(aes(ymin=1, ymax = top, xmin=0, xmax=Ti, fill="during ART"), alpha = .8) +
    geom_rect(aes(ymin=1, ymax = top, xmin=Ti, xmax = t_end, fill="after ART stop"), alpha = .8) +
    scale_fill_brewer(palette = 'Pastel2', name = 'Treatment Stage') +
    theme_classic() +
    # plot data
    geom_point(data = subdat, aes(x = time_years, y = log_rna_halfDL, colour = "log10(RNA)", shape = factor(censor))) +
    geom_point(data = subdat_cd4, aes(x = time_CD4, y = y_cd4_recale, colour = "log(CD4)")) +
    geom_line(data = subdat, aes(x = time_years, y = log_rna_halfDL, colour = "log10(RNA)")) +
    geom_line(data = subdat_cd4, aes(x = time_CD4, y = y_cd4_recale, colour = "log(CD4)")) +
    geom_hline(yintercept = log10(40), linetype = 'dashed') +
    geom_text(data = subdat, x = 9/10*t_end, y = 1.4, label = "detection limit") +
    geom_vline(xintercept = upper_bound, linetype = 'dotted') +
    geom_text(data = subdat, x = 1.05*upper_bound, y = 1.1*y, label = "Latest\nchange\npoint") +
    geom_vline(xintercept = t_max_cd4, linetype = 'dotted') +
    geom_text(data = subdat, x = 0.95*t_max_cd4, y = 0.8*y, label = "Possible\nchange\npoint") +
    scale_shape_manual(name = "Censor", values=c(19, 1), labels = c("No", "Yes")) +
    scale_colour_manual(name = "Biomarker",
                        values = c("blue", "black")) +
    xlab("Time in years") +
    scale_y_continuous(
      
      # Features of the first axis
      name = bquote("Viral load (in" ~ log[10]~"-scale)"),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(~ rescale(., from = c(min(subdat$log_rna_halfDL), max(subdat$log_rna_halfDL)), 
                                    to = c(min(subdat_cd4$y_logCD4), max(subdat_cd4$y_logCD4))), 
                          name=bquote("CD4 count (in" ~ log[10]~"-scale)"))
    )
  p
  ggsave(paste0("Ind_plot_cleaned_paper/ID", i, ".png"), width = 10, height = 6)
  cat("i=", i, '\n')
}


i = 32
subdat <- data %>% filter(PATIENT == i)
subdat_cd4 <- dplyr::filter(data_cd4, PATIENT == i)
subdat_cd4$y_cd4_recale <- rescale(subdat_cd4$y_logCD4, 
                                   from = c(min(subdat_cd4$y_logCD4), max(subdat_cd4$y_logCD4)), 
                                   to = c(min(subdat$log_rna_halfDL), max(subdat$log_rna_halfDL)))

t_end <- max(subdat$time_years)
Ti <- unique(subdat$treatment_stop_years)
upper_bound <- unique(subdat$first_reboundobs)
subdat_2 <- subdat %>% filter(censor == 0)
y <- median(subdat_2$log_rna)
top <- 1.1*max(subdat$log_rna_halfDL)
treatment_stop <- unique(subdat$treatment_stop_years)
subdat_cd4_2 <- subdat_cd4 %>% filter(time_CD4 <= upper_bound & time_CD4 >= treatment_stop)
t_max_cd4 <- max(subdat_cd4_2[which(subdat_cd4_2$y_logCD4 == max(subdat_cd4_2$y_logCD4)), ]$time_CD4)
subdat$measure <- "Viral load"
subdat_cd4$measure <- "CD4"
p <- ggplot() +
  xlim(0, t_end) +
  ylim(0, top) + 
  theme_classic() +
  theme(text = element_text(size = 16))+
  
  # Viral load: circle (19) and open circle (1 for censored)
  geom_point(data = subdat,
             aes(x = time_years, y = log_rna_halfDL,
                 shape = factor(censor)), size = 3) +
  
  # CD4: triangle (17)
  geom_point(data = subdat_cd4,
             aes(x = time_CD4, y = y_cd4_recale,
                 shape = measure), size = 3) +
  
  geom_line(data = subdat, aes(x = time_years, y = log_rna_halfDL)) +
  geom_line(data = subdat_cd4, aes(x = time_CD4, y = y_cd4_recale)) +
  
  geom_hline(yintercept = log10(40), linetype = 'dashed') +
  geom_text(data = subdat, x = 9/10*t_end, y = 1.4, size = 6, label = "detection limit") +
  geom_vline(xintercept = upper_bound, linetype = 'dotted') +
  geom_text(data = subdat, x = upper_bound, y = 2.5, size = 6, label = "Latest\nchange\npoint") +
  geom_vline(xintercept = Ti, linetype = 'dotted') +
  geom_text(data = subdat, x = Ti, y = 2.5, size = 6, label = "Treatment\nstop\ntime") +
  geom_vline(xintercept = t_max_cd4, linetype = 'dotted') +
  geom_text(data = subdat, x = t_max_cd4, y =3.2, size = 6, label = "Possible\nchange\npoint") +
  
  # manual shape control:
  scale_shape_manual(
    name = "Measurement",
    values = c("0" = 19, "1" = 1, "CD4" = 17),  # 17 = triangle
    labels = c("0" = "Viral load", "1" = "Viral load (censored)", "CD4" = "CD4")
  ) +
  
  xlab("Time in years") +
  scale_y_continuous(
    name = bquote("Viral load (in" ~ log[10]~"-scale)"),
    sec.axis = sec_axis(~ rescale(., from = c(min(subdat$log_rna_halfDL),
                                              max(subdat$log_rna_halfDL)),
                                  to   = c(min(subdat_cd4$y_logCD4),
                                           max(subdat_cd4$y_logCD4))),
                        name=bquote("CD4 count (in" ~ log[10]~"-scale)"))
  )
p
ggsave(paste0("Ind_plot_cleaned_paper/ID", i, "_v2.png"), width = 20, height = 6)
cat("i=", i, '\n')

data_decay <- data %>% 
  dplyr::select(PATIENT, time_years, log_rna, log_rna_halfDL, censor, treatment_stop_years, treatment, last_decayobs, first_reboundobs)
colnames(data_decay) <- c("PATIENT", "time_decay", "y_decay", "y_decay_halfDL", "censor_decay", "treatment_stop", "treatment", "last_decayobs", "first_reboundobs")

data_rebound <- data %>% 
  dplyr::select(PATIENT, time_years, log_rna, log_rna_halfDL, censor, treatment_stop_years, treatment, last_decayobs, first_reboundobs)
colnames(data_rebound) <- c("PATIENT", "time_rebound", "y_rebound", "y_rebound_halfDL", "censor_rebound", "treatment_stop", "treatment", "last_decayobs", "first_reboundobs")

data_cd4 <- data_cd4 %>% 
  dplyr::select(PATIENT, time_CD4, y_CD4, y_logCD4)

data_trans <- data_rebound %>% 
  filter(treatment == 0) %>% 
  dplyr::select(PATIENT, last_decayobs, first_reboundobs, treatment_stop) %>% 
  group_by(PATIENT) %>%
  slice(1) %>%
  ungroup()

write.csv(data_decay, here::here("Data", "cleaned_data_decay.csv"), row.names = FALSE)
write.csv(data_rebound, here::here("Data", "cleaned_data_rebound.csv"), row.names = FALSE)
write.csv(data_cd4, here::here("Data", "cleaned_data_cd4.csv"), row.names = FALSE)
write.csv(data_trans, here::here("Data", "cleaned_data_trans.csv"), row.names = FALSE)

