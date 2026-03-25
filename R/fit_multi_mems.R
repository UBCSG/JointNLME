library(matrixcalc)
library(Deriv)
library(plyr)

#' Fit multiple mixed-effect models jointly based on h-likelihood method.
#'
#' This function is for joint modelling of multiple mixed-effect models.
#' The continuous longitudinal data may be left-censored due to lower
#' limits of quantification.
#' @import matrixcalc
#' @import magrittr
#' @importFrom plyr compact
#' @useDynLib hhjm, .registration = TRUE
#' @param pre_memList a list of sublists. Each sublist provides necessary
#' information of an LME, GLME, or NLME model for pre-data;
#' @param post_memList a list of sublists. Each sublist provides necessary
#' information of an LME, GLME, or NLME model for post-data;
#' @param pre_data pre ART-interruption data;
#' @param post_data post ART-interruption data;
#' @param subject_id character; name of subject id;
#' @param randeff_info assumptions of random effects; distribution = c("normal", "t-dist");
#' degree is used only when distribution = "t-dist"; correlation = NULL or "independent".
#' @param loglike_tol convergence criteria based on changes in log-likelihood values;
#' @param par_tol convergence criteria based on changes in parameter estimates;
#' @param iterMax maximum number of iterations;
#' @param Silent TRUE or FALSE; If FALSE, print out messages;
#' @export
fit_multi_mems <- function(memList, dataList, subject_id,
                           randeff_info = list(distribution = "normal", degree = 0, correlation = NULL),
                           cov_method = c("spherical", "cholesky"),
                           loglike_tol = 1e-3, par_tol = 1e-2, iterMax = 10,
                           REML = FALSE, Silent = TRUE, naive = FALSE, adjust = FALSE){

  check_memList(memList)
  cov_method <- match.arg(cov_method)
  RespLog <- lapply(memList, function(x){get_log_density(x)$log_density})
  random_effects <- extract_random_effects(memList)
  fixed_param_info <- extract_fixed_param(memList)
  disp_param_info <- extract_disp_param(memList)
  invSIGMA <- initiate_invSIGMA(memList, randeff_info)
  
  fixed_param <- fixed_param_info$fixed_param
  disp_param <- disp_param_info$disp_param
  uniqueID <- lapply(dataList, function(dat){dat[, subject_id]}) %>% unlist %>% unique()
  fix_pars <- fixeff_par(RespLog, random_effects, fixed_param, memList)
  disp_pars <- disp_par(RespLog, c(random_effects, names(fixed_param)), disp_param)
  likDiff <- Diff <- convergence <- m <- 1
  trans_index <- which(sapply(memList, function(x) x$event == "trans") == TRUE)

  while(likDiff > loglike_tol & Diff > par_tol & m < iterMax){
    ################################################
    ########### Estimate random effects ############
    ################################################
    # estimate random effects by max log h-likelihood
    cat("############## Iteration:", m, "###############","\n")
    cat("# Estimating random effects")
    

    # Given dataList which should contain estimates of baseline hazard  but not random effects
    BiTi_df_list <- est_randeff(uniqueID, RespLog, dataList, subject_id,
                                random_effects, invSIGMA,
                                par_val = c(disp_param, fixed_param),
                                distribution = randeff_info$distribution,
                                degree = randeff_info$degree,
                                Silent = Silent, Scale = (cov_method == "spherical"), 
                                mc.cores = 1, naive)
    
    Bi_df <- BiTi_df_list$Bi_df_center
    Ti_df <- BiTi_df_list$Ti_df
    
    # new_dataList contains estimates of baselinehazard and random effects
    new_dataList <- lapply(dataList, function(dat){
      dat <- dplyr::left_join(dplyr::select(dat, - tidyselect::one_of(random_effects, "Ti0_sim")), Bi_df, by = subject_id)
      dplyr::left_join(dat, Ti_df, by = subject_id)
    })
    
    new_dataList[[1]] <- new_dataList[[1]] %>%
      filter(time_decay < Ti0_sim + treatment_stop)
    new_dataList[[2]] <- new_dataList[[2]] %>%
      filter(time_rebound >= Ti0_sim + treatment_stop) %>%
      mutate(time_rebound = time_rebound - (Ti0_sim + treatment_stop))
    
    new_dataList[[trans_index]] <- new_dataList[[trans_index]] %>%
      mutate(time = Ti0_sim) %>% 
      mutate(time = ifelse(time == last_decayobs, last_decayobs + 0.01, time)) %>% 
      mutate(time = ifelse(time == first_reboundobs, first_reboundobs - 0.01, time))
    
    new_dataList[[trans_index]] %>% filter(time < last_decayobs | time > first_reboundobs)
    
    cat('\n')

    ################################################
    ############# Estimate parameters ##############
    ################################################
    # Non-robust estimate of fixed parameters
    cat("# Estimating fixed parameters")
    new_fixed_param <- est_fixeff(param = fixed_param,
                                 other_param = disp_param,
                                 lower = fixed_param_info$lower,
                                 upper = fixed_param_info$upper,
                                  RespLog, new_dataList,
                                  Bi_df, invSIGMA, 
                                  distribution = randeff_info$distribution,
                                  degree = randeff_info$degree,
                                  fix_pars = fix_pars, subject_id,
                                  all_glmeObjects = memList,
                                  Silent = Silent, naive, adjust)
    
    cat("\n")
    
    # # Non-robust estimate of dispersion parameters
    cat("# Estimating dispersion parameters ...")
    res_disp <- est_disp(param = disp_param,
                             other_param = new_fixed_param,
                             RespLog,
                             new_dataList, Bi_df, invSIGMA = invSIGMA, 
                             lower = disp_param_info$lower,
                             upper = disp_param_info$upper,
                             distribution = randeff_info$distribution,
                             degree = randeff_info$degree,
                             correlation = randeff_info$correlation,
                             cov_method = cov_method,
                             fix_pars = fix_pars,
                             REML = REML,
                             all_glmeObjects = memList,
                             Silent = Silent, naive, adjust)
    cat("\n")
    new_disp_param <- res_disp$disp_param
    new_invSIGMA <- res_disp$invSIGMA

    
    ####################################################
    ################## update results ##################
    ####################################################
    # calculate approximated log marginal likelihood value
    
    if (adjust == TRUE){
      new_loglike_value <- approx_log_like(RespLog, par_val = c(new_fixed_param, new_disp_param),
                                           new_dataList, Bi_df,
                                           invSIGMA = new_invSIGMA, 
                                           distribution = randeff_info$distribution,
                                           degree = randeff_info$degree,
                                           all_glmeObjects = memList,
                                           fix_pars = fix_pars, naive)
    } else {
      new_loglike_value <- eval_log_hlike(RespLog, new_dataList, par_val = c(new_fixed_param, new_disp_param), 
                                          Bi = as.matrix(Bi_df[, -1]), invSIGMA = new_invSIGMA, 
                                          distribution = randeff_info$distribution, 
                                          degree = randeff_info$degree)
    }

    if(m == 1){
      likDiff <- 1
    } else{
      likDiff <- abs((new_loglike_value - loglike_value)/loglike_value)
    }

    # calcuate relative changes in mean parameters
    old_fixed_param <- c(fixed_param, disp_param)[names(new_fixed_param)]
    old_disp_param <- c(fixed_param, disp_param)[names(new_disp_param)]
    relative_change <- abs((unlist(new_fixed_param) - unlist(old_fixed_param))/(unlist(old_fixed_param) + 1e-16))
    Diff <- mean(relative_change)

    # print iterating result
    print("Approximate log likelihood:")
    print(new_loglike_value)
    print("Estimates of fixed parameters:")
    print(round(unlist(new_fixed_param), 3))
    print("Estimates of dispersion parameters:")
    print(round(unlist(new_disp_param), 3))
    print("Maximum relative changes in fixed parameters:")
    print(round(Diff, 4))
    print("Relative change in log likelihood:")
    print(likDiff)
    cat("##########################################","\n")

    fixed_param <- new_fixed_param
    invSIGMA <- new_invSIGMA
    disp_param <- new_disp_param
    loglike_value <- new_loglike_value
    m <- m + 1
  } 

  # messages about convergence success or failure
  if((likDiff > loglike_tol  & Diff > par_tol)){
    warning("Iteration limit reached without covergence.")
    convergence <- 1
  }
  if(likDiff <= loglike_tol & likDiff >= 0){
    message("Successful convergence. Iteration stops because likDiff <= loglike_tol.")
    convergence <- 0
  }
  if(Diff <= par_tol){
    message("Successful convergence. Iteration stops because FixedParDiff <= par_tol.")
    convergence <- 0
  }

  fixed_param_sd <- NULL
  if(convergence == 0){
    fixed_param_sd <- get_sd_hlike(fixed_param, random_effects,  other_param = disp_param,
                             RespLog, new_dataList, Bi_df, invSIGMA, 
                             disp_pars, 
                             all_glmeObjects = memList,
                             distribution = randeff_info$distribution,
                             degree = randeff_info$degree, naive)
    
    print("Estimates of fixed parameters SD:")
    print(round(unlist(fixed_param_sd), 4))
  }

  output <- list(fixed_est = unlist(fixed_param),
                 fixed_sd = fixed_param_sd,
                 disp_est = unlist(disp_param),
                 Bi = Bi_df,
                 Ti = Ti_df,
                 invSIGMA = invSIGMA,
                 random_effects = random_effects,
                 convergence = convergence,
                 loglike_value = loglike_value,
                 RespLog = RespLog,
                 dataList = dataList,
                 memList = memList,
                 subject_id = subject_id,
                 randeff_info = randeff_info)
  class(output) <- "hhjm"
  output
}


