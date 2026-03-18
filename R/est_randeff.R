
#' Estimate random effects by maximizing h-likelihood function,
#' given fixed and dispersion parameters
#'
est_randeff_by_hlike <- function(RespLog, sub_dataList, 
                                 random_effects, invSIGMA, par_val, 
                                 distribution, degree, Silent = T, naive){
  q <- length(random_effects)
  
  ff <- function(xx, ...){
    fy <- numeric(1)
    
    BiTi <- data.frame(Vassign(c(random_effects, "Ti0_sim"), xx))
    Bi <- BiTi[, 1:q]
    
    new_sub_dataList <- lapply(sub_dataList, function(dat){merge(dat, BiTi)})
    new_sub_dataList[[1]] <- new_sub_dataList[[1]] %>%
      mutate(Ti_sim = Ti0_sim + treatment_stop) %>%
      filter(time_decay < Ti_sim)
    new_sub_dataList[[2]] <- new_sub_dataList[[2]] %>%
      mutate(Ti_sim = Ti0_sim + treatment_stop) %>%
      filter(time_rebound >= Ti_sim) %>%
      mutate(time_rebound = time_rebound - Ti_sim)
    if (naive == FALSE){
      new_sub_dataList[[3]]$fitted_CD4 <- eval_fn_row(Object_CD4$reg_equation, new_sub_dataList[[3]], par_val, get_gradient = F)
      
      if (nrow(new_sub_dataList[[4]]) != 0){
        new_sub_dataList[[4]] <- new_sub_dataList[[4]] %>%
          mutate(time = Ti0_sim) %>%
          mutate(maxCD4 = max(new_sub_dataList[[3]]$fitted_CD4))
      }
    } else {
      if (nrow(new_sub_dataList[[3]]) != 0){
        new_sub_dataList[[3]] <- new_sub_dataList[[3]] %>%
          mutate(time = Ti0_sim)
      }
    }
    
    
    fy <- eval_log_hlike(RespLog, new_sub_dataList, par_val, as.matrix(Bi),
                         invSIGMA, distribution, degree)
    
    - fy
  }
  
  str_bi <- rmvnorm(1, rep(0, q), solve(invSIGMA))
  str_Bi <- data.frame(Vassign(c(random_effects), str_bi))
  
  if (naive == FALSE){
    sub_dataList[[3]]$fitted_CD4 <- eval_fn_row(Object_CD4$reg_equation, sub_dataList[[3]], c(par_val, str_Bi), get_gradient = F)
    # str_Ti <- exp(par_val$omega1 * str_Bi$tau2i + par_val$omega2 * max(sub_dataList[[3]]$fitted_CD4))
    lower_new <- c(rep(-Inf, q), sub_dataList[[4]]$last_decayobs)
    upper_new <- c(rep(Inf, q), sub_dataList[[4]]$first_reboundobs - sub_dataList[[4]]$treatment_stop)
    str_Ti <- upper_new[q + 1] * 0.5
  } else{
    # str_Ti <- exp(par_val$omega1 * str_Bi$tau2i + par_val$omega2 * sub_dataList[[3]]$maxCD4_obs)
    lower_new <- c(rep(-Inf, q), sub_dataList[[3]]$last_decayobs)
    upper_new <- c(rep(Inf, q), sub_dataList[[3]]$first_reboundobs - sub_dataList[[3]]$treatment_stop)
    str_Ti <- upper_new[q + 1] * 0.5
  }
   
  res <- optim_iterate_fix(
    str_val = c(str_bi, str_Ti),
    ff, gr = NULL, 
    check = 1 - as.numeric(Silent),
    gr_tol = 1e-3, random_start = TRUE,
    lower = lower_new,
    upper = upper_new,
    # extra info
    RespLog = RespLog, 
    # gr_logh_b = gr_logh_b,
    sub_dataList = sub_dataList,
    random_effects = random_effects,
    invSIGMA = invSIGMA, par_val = par_val,
    distribution = distribution, degree = degree,
    q = q,
    e4 = e4)
  res$gamma
}

#' Obtain estimates of random effects
#'
#' @importFrom dplyr select
#' @importFrom tidyselect one_of
est_randeff <- function(uniqueID,
                        RespLog, dataList, subject_id, 
                        random_effects, invSIGMA, par_val, 
                        distribution = "normal", degree = 0,
                        Silent = TRUE, Scale = FALSE, mc.cores = 1, naive){
  q <- length(random_effects)
  stopifnot(length(RespLog) == length(dataList))
  if(class(par_val) != "list"){
    par_val <- as.list(par_val)
  }
  
  BiTi_est <- parallel::mclapply(uniqueID, function(i){
    cat(".")
    
    sub_dataList <- lapply(dataList, function(dat){
      subset(dat, dat[, subject_id] == i) %>% dplyr::select(-tidyselect::one_of(c(random_effects, "Ti0_sim", "maxCD4", "time")))
    })
    
    biTi <- try(est_randeff_by_hlike(RespLog, sub_dataList, random_effects,
                                     invSIGMA, par_val, distribution, degree, Silent = T, naive), silent = T)
    if(class(biTi) == "try-error"){
      stop(paste("Estimating random effects failed for subject", i))
    }
    biTi
  }, mc.cores = mc.cores) %>% do.call(rbind, .)
  
  Bi_est = BiTi_est[, 1:q]
  
  Bi_center <- scale(Bi_est, center = TRUE, scale = Scale)
  if (length(uniqueID) == 1){
    Bi_df_center <- cbind(uniqueID, data.frame(t(Bi_center)))
    Bi_df <- cbind(uniqueID, data.frame(t(Bi_est)))
  } else {
    Bi_df_center <- cbind(uniqueID, data.frame(Bi_center))
    Bi_df <- cbind(uniqueID, data.frame(Bi_est))
  }
  names(Bi_df) <- names(Bi_df_center) <- c(subject_id, random_effects)
  
  Ti_df <- cbind(uniqueID, data.frame(BiTi_est[, q+1]))
  names(Ti_df) <- c(subject_id, "Ti0_sim")
  
  list(Bi_df = Bi_df, Bi_df_center = Bi_df_center, Ti_df = Ti_df)
  # list(Bi_df = Bi_df, Bi_df_center = Bi_df_center)
}
