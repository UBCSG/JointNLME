library(ggplot2)
library(magrittr)
library(lme4)
library(readxl)
library(saemix)
library(dplyr)
library(nlme)
library(Deriv)
library(matrixcalc)
library(optimx)

setwd(here::here("R"))
file.sources = list.files(pattern = "*.R")
sapply(file.sources, source, .GlobalEnv)


setwd(here::here("Data"))
### input data sets
pre_art_data_cd
4 <- read.csv("ex_pre_art_data_cd4.csv")
pre_art_data_rna <- read.csv("ex_pre_art_data_rna.csv")
post_art_data <- read.csv("ex_post_art_data.csv")

### detection limit
limit_value = log10(40)

## Form models
glmeObject_CD4 <- list(
  modeltype = "nlme",
  response = "log_cd4",
  reg_equation = "alpha0 + d11*b11 + (alpha1 + d12*b12) * time_k",
  distribution = 'normal',
  random_effects = c('b11', "b12"),
  fixed_param = list(names = paste0("alpha", 0:1),
                     start_values = c(6.2, 0.4)),
  disp_param = list(names = c("d11", "d12"),
                    start_values = c(1, 0.1),
                    lower = c(0, 0),
                    upper = c(Inf, Inf)),
  sigma = 0.5)

pre_glmeObject <- list(
  modeltype = "nlme",
  response = "rna",
  reg_equation = "log10(exp(eta0 + d21 * b21) + exp(eta1 + d22*b22 - exp(eta2 + d23*b23)*time_m))",
  distribution = 'normal',
  random_effects = c("b21", "b22", "b23"),
  fixed_param = list(names = paste0("eta", 0:2),
                     start_values = c(2, 11, 1)),
  disp_param = list(names = c("d21", "d22", "d23"),
                    start_values = c(0.1, 0.1, 0.1),
                    lower = c(0, 0, 0),
                    upper = c(Inf, Inf, Inf)),
  sigma = 0.5,
  left_censoring = list(indicator = "z",
                        limit_value = limit_value,
                        method = "tobit"))


post_glmeObject <- list(
  modeltype = "nlme",
  response = "log_rna",
  reg_equation = "(beta1 + d31*b31)* time / (time + exp(beta2 - (beta3)*time)) + (beta4)",
  distribution = 'normal',
  random_effects = c("b31"),
  fixed_param = list(names = c("beta1", "beta2", "beta3", "beta4"),
                     start_values = c(3, 8, 3, 1)),
  disp_param = list(names = c("d31"),
                    start_values = c(0.5),
                    lower = c(0),
                    upper = c(Inf)),
  sigma = 0.5,
  left_censoring = list(indicator = "c",
                        limit_value = limit_value,
                        method = "tobit"))

tran_glmeObject <- list(
  modeltype = "survival",
  event = "rebound",
  reg_equation = "(gamma1 * b21 + gamma2 * b22 + gamma3 * b23 + gamma4 * b11 + gamma5 * b12)",
  distribution = NULL,
  random_effects = NULL,
  fixed_param = list(names = c("gamma1", "gamma2", "gamma3", "gamma4","gamma5"),
                     start_values = c(3, 8, 3, 1)),
  disp_param = NULL,
  sigma = 0.5)

memList <- list(glmeObject_CD4, pre_glmeObject, post_glmeObject, tran_glmeObject)
dataList <- list(pre_art_data_cd4, pre_art_data_rna, post_art_data)

md0 <- fit_multi_mems(memList = memList,
                      dataList = dataList,
                      subject_id = "PATIENT",
                      loglike_tol = 1e-2, par_tol = 1e-2,
                      Silent = T, iterMax = 30)

summary.hhjm(md0)
