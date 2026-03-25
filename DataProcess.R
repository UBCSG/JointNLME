# =============================================================================
# HIV/AIDS Joint Modeling: Data Processing
# =============================================================================

library(tidyverse)
library(ggpubr)
library(berryFunctions)

setwd(here::here())

# Source all helper functions from the R/ directory
file.sources <- list.files("R", pattern = "\\.R$", full.names = TRUE)
sapply(file.sources, source, .GlobalEnv)

# =============================================================================
# 1. DATA IMPORT & CLEANING
# =============================================================================

raw_data <- read.csv(here::here("Data", "raw_data.csv"), header = TRUE)

# Remove observations missing seroconversion date or viral RNA
data_all <- raw_data %>%
  filter(!is.na(days_from_seroco), !is.na(RNA_V)) %>%
  mutate(
    log_cd4           = log(CD4_V),                                  # NA preserved automatically
    censor            = as.integer(RNA_V <= 40),                     # 1 = below detection limit
    log_rna           = ifelse(censor == 1, log10(40), log10(RNA_V)),
    log_rna_halfDL    = ifelse(RNA_V <= 40, log10(20), log10(RNA_V)),# half detection-limit substitution
    time_days         = days_from_seroco - art_start,
    time_years        = time_days / 365,
    treatment_stop    = ati_start - art_start,
    treatment_stop_years = treatment_stop / 365,
    limit             = 0
  ) %>%
  filter(days_from_seroco >= art_start) %>%   # keep only on-ART observations
  arrange(PATIENT, time_days)

# Attach each patient's peak CD4 value (used downstream)
max_cd4 <- data_all %>%
  filter(!is.na(log_cd4)) %>%
  group_by(PATIENT) %>%
  slice_max(log_cd4, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(PATIENT, time_log_cd4_max = days_from_seroco, log_cd4_max = log_cd4)

data_all <- data_all %>%
  left_join(max_cd4, by = "PATIENT") %>%
  mutate(last_decayobs = 0)


# =============================================================================
# 2. MANUAL OUTLIER REMOVAL
# =============================================================================
# Observations removed based on clinical review (see Methods section).
# Each row removes either a single time-point or all data beyond a threshold.

outliers_exact <- tribble(
  ~PATIENT, ~time_days,
  3,  875,  5, 1040,   7,  389,
  11, 670,  11,  706,  12,  732,  12,  942,
  14, 784,  16,  630,  16,  679,  17,  698,
  21, 550,  21,  579,  21,  818,  21,  865,
  24, 454,  29,  658,  32,  503,  32,  559,
  46, 749,  53,  445,  56,  445,  57,  982,
  73, 210,  73,  232,  75,  623,  76,  786
)

outliers_threshold <- tribble(
  ~PATIENT, ~max_days,
  2,  1400,  3,  1100,  5,  1100,  7,   900,
  8,  1800,  9,   900, 11,  1200, 12,  1200,
  14,  900, 15,  1100, 16,  1200, 17,  1000,
  18, 1500, 19,   800, 20,  1200, 21,   750,
  23,  450, 24,  1000, 25,  1000,
  28, 1100, 29,  1000, 30,  1200,
  33, 1000, 35,  1000, 36,   800, 37,  1000,
  38, 1200, 40,  1200, 42,  1000, 43,  1200,
  47,  800, 50,  1500, 52,   900, 53,  1000,
  56, 1000, 57,  1500, 61,  1100, 64,  1000,
  66, 1200, 67,   800, 69,   800, 70,   800,
  73, 1000, 75,  1200, 76,  1200
)

# Patients excluded entirely (no usable data after review)
patients_excluded <- c(4, 26, 27, 41, 62, 71)

data <- data_all %>%
  anti_join(outliers_exact, by = c("PATIENT", "time_days")) %>%
  anti_join(outliers_threshold %>%
              dplyr::rename(threshold = max_days) %>%
              inner_join(data %>% dplyr::select(PATIENT, time_days), by = "PATIENT") %>%
              filter(time_days > threshold) %>%
              dplyr::select(PATIENT, time_days),
            by = c("PATIENT", "time_days")) %>%
  filter(!PATIENT %in% patients_excluded)


# =============================================================================
# 3. FIRST REBOUND OBSERVATION
# =============================================================================
# For each patient, identify the earliest detectable viral load after ART stop.
# Some values are overridden manually based on clinical judgment.

first_reboundobs <- data %>%
  filter(treatment == 0, !censor) %>%
  group_by(PATIENT) %>%
  slice_min(time_years, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(PATIENT, first_reboundobs = time_years)

manual_rebound_overrides <- tribble(
  ~PATIENT, ~first_reboundobs,
  11,  782 / 365,
  18,  752 / 365,
  26,  348 / 365,
  32, 1042 / 365,
  41,  296 / 365,
  56,  581 / 365,
  62,  327 / 365,
  71,  624 / 365
)

data <- data %>%
  left_join(first_reboundobs, by = "PATIENT") %>%
  rows_update(manual_rebound_overrides, by = "PATIENT", unmatched = "ignore") %>%
  arrange(PATIENT, time_years)

write.csv(data, here::here("Data", "cleaned_data.csv"), row.names = FALSE)

# =============================================================================
# 4. EXPORT ANALYSIS-READY DATASETS
# =============================================================================

common_cols <- c("PATIENT", "time_years", "log_rna", "log_rna_halfDL",
                 "censor", "treatment_stop_years", "treatment",
                 "last_decayobs", "first_reboundobs")

data_decay <- data %>%
  dplyr::select(all_of(common_cols)) %>%
  dplyr::rename(time_decay    = time_years,
         y_decay       = log_rna,
         y_decay_halfDL = log_rna_halfDL,
         censor_decay  = censor,
         treatment_stop = treatment_stop_years)

data_rebound <- data %>%
  dplyr::select(all_of(common_cols)) %>%
  dplyr::rename(time_rebound     = time_years,
         y_rebound        = log_rna,
         y_rebound_halfDL = log_rna_halfDL,
         censor_rebound   = censor,
         treatment_stop   = treatment_stop_years)

data_cd4 <- data_all %>%
  dplyr::select(PATIENT, time_years, CD4_V, log_cd4) %>% 
  filter(!is.na(log_cd4))
colnames(data_cd4)[2:4] <- c("time_CD4", "y_CD4", "y_logCD4")

# One row per patient: transition-window summary for the rebound model
data_trans <- data_rebound %>%
  filter(treatment == 0) %>%
  group_by(PATIENT) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(PATIENT, last_decayobs, first_reboundobs, treatment_stop)

write.csv(data_decay,   here::here("Data", "cleaned_data_decay.csv"),   row.names = FALSE)
write.csv(data_rebound, here::here("Data", "cleaned_data_rebound.csv"), row.names = FALSE)
write.csv(data_cd4, here::here("Data", "cleaned_data_cd4.csv"),     row.names = FALSE)
write.csv(data_trans,   here::here("Data", "cleaned_data_trans.csv"),   row.names = FALSE)
