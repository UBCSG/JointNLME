#' Estimate Ti by maximizing the profile h-likelihood function
#'
#' @importFrom dplyr select left_join
#' @importFrom Deriv Simplify
#' @importFrom matrixcalc matrix.trace
est_Ti <- function(RespLog, sub_dataList, 
                   Bi_df, invSIGMA, par_val, 
                   distribution, degree, last_decayobs, Silent = T){
  
  Ti_df <- Bi_df %>% 
    mutate(logTi = par_val$omega1 * tau2i + par_val$omega2 * a1i + rnorm(1, 0, par_val$time_sigma)) %>% 
    dplyr::select(ID, logTi)
  
  exp_dataList <- lapply(new_dataList, function(dat){
    dplyr::left_join(dplyr::select(dat, - tidyselect::one_of("logTi")), Ti_df, by = subject_id)
  })
  
  # Update the data with the simulated time point
  exp_dataList[[1]] <- exp_dataList[[1]] %>%
    mutate(Ti_sim = exp(logTi)) %>%
    filter(time_decay <= Ti_sim + last_decayobs)
  exp_dataList[[2]] <- exp_dataList[[2]] %>%
    mutate(Ti_sim = exp(logTi)) %>%
    filter(time_rebound > Ti_sim + last_decayobs) %>%
    mutate(time_rebound = time_rebound - (Ti_sim + last_decayobs))
  exp_dataList[[3]] <- exp_dataList[[3]] %>%
    mutate(Ti_sim = exp(logTi)) %>%
    filter(time <= Ti_sim + last_decayobs)
  exp_dataList[[4]] <- exp_dataList[[4]] %>%
    mutate(Ti_sim = exp(logTi)) %>%
    mutate(trans = ifelse(time <= (Ti_sim + last_decayobs), 0, 1)) %>%
    group_by(ID) %>%
    slice(min(which(trans == 1))) %>%
    as.data.frame() %>%
    mutate(time = Ti_sim + last_decayobs)
  
  
  logTi <- data.frame(Vassign("logTi", xx))
  
  # Negative adjusted profile h-likelihood with param as input
  ff <- function(xx, ...){
    # all parameter values
    par_val <- Vassign(names(param), xx) %>% c(., other_param)
    
    lhlike <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]), invSIGMA, distribution, degree)
    - lhlike
  }
  
  lapply(sub_dataList, function(dat){
    new_sub_dataList <- sub_dataList %>% 
      mutate(e4 = rnorm(1, 0, par_val$time_sigma)) %>% 
      mutate(logTi = par_val$omega1 * tau2i + par_val$omega2 * a1i + e4)
  })
  
 
  
  
  ## gr() function is not used because it's slower than calculating numerical gradients
  res <- optim_iterate_fix(str_val = unlist(param), ff, gr = NULL, lower = lower, upper = upper,
                           check = 1 - as.numeric(Silent),
                           # methods = "nlminb",
                           # extra inputs
                           param = param, other_param = other_param,
                           RespLog = RespLog,
                           Bi_df = Bi_df, invSIGMA = invSIGMA,
                           distribution = distribution, degree = degree,
                           random_effects = random_effects,
                           new_dataList = new_dataList,
                           Hmats = Hmats,
                           negH3 = negH3,
                           dH_theta = dH_theta,
                           dfs = dfs,
                           uniqueID = uniqueID, fix_pars = fix_pars,
                           all_glmeObjects = all_glmeObjects)
  res$gamma
}


