#' Estimate fixed parameters by maximizing the profile h-likelihood function
#'
#' @importFrom dplyr select left_join
#' @importFrom Deriv Simplify
#' @importFrom matrixcalc matrix.trace
est_fixeff_sim <- function(param, other_param,
                           lower = -Inf, upper = Inf,
                           RespLog, new_dataList, 
                           Bi_df, invSIGMA, 
                           distribution = "normal", degree = 0, 
                           fix_pars = fix_pars, subject_id,
                           all_glmeObjects, 
                           Silent = TRUE, naive, adjust){
  subject_id <- names(Bi_df)[1]
  random_effects <- names(Bi_df)[-1]
  q <- length(random_effects)
  n <- nrow(Bi_df)
  p <- length(param)
  uniqueID <- unique(Bi_df[, 1])
  Hmats <- fix_pars$Hmats
  
  if(distribution == "normal"){
    negH3 <- lapply(1:n, function(i){ -invSIGMA })
  } else if (distribution == "t-dist"){
    numo <- as.matrix(Bi_df[, -1]) %*% invSIGMA
    deno <- numo %*% t(Bi_df[, -1]) %>% diag()
    negH3 <- lapply(1:n, function(i){
      numo1 <- (degree + q)*invSIGMA*(degree + deno[i]) -
        (degree + q)*as.matrix(numo[i, ])%*%(2*numo[i, ])
      -numo1/((degree + deno[i])^2)
    })
  }
  
  # Negative adjusted profile h-likelihood with param as input
  ff_adj <- function(xx, ...){
    # all parameter values
    par_val <- Vassign(names(param), xx) %>% c(., other_param)
    
    if (naive == FALSE){
      new_dataList[[3]]$fitted_CD4 <- eval_fn_row(Object_CD4$reg_equation, new_dataList[[3]], par_val, get_gradient = F)
      max_CD4 <- new_dataList[[3]] %>%
        group_by(!!sym(subject_id)) %>%
        filter(fitted_CD4 == max(fitted_CD4)) %>%
        dplyr::select(!!sym(subject_id), fitted_CD4) 
      colnames(max_CD4)[2] <- "maxCD4"
      new_dataList[[4]] <- new_dataList[[4]] %>% dplyr::select(-tidyselect::one_of(c("maxCD4")))
      new_dataList[[4]] <- merge(new_dataList[[4]], max_CD4) 
      exp_dataList <- new_dataList
      for(j in 1:length(all_glmeObjects)){
        exp_dataList[[j]][all_glmeObjects[[j]]$response] <-
          eval_fn_row(all_glmeObjects[[j]]$reg_equation, exp_dataList[[j]], par_val, get_gradient = F)
      }
      exp_dataList[[4]][all_glmeObjects[[4]]$response] <- exp(exp_dataList[[4]][all_glmeObjects[[4]]$response])
    } else{
      exp_dataList <- new_dataList
      for(j in 1:length(all_glmeObjects)){
        exp_dataList[[j]][all_glmeObjects[[j]]$response] <-
          eval_fn_row(all_glmeObjects[[j]]$reg_equation, exp_dataList[[j]], par_val, get_gradient = F)
      }
      exp_dataList[[3]][all_glmeObjects[[3]]$response] <- exp(exp_dataList[[3]][all_glmeObjects[[3]]$response])
    }
    
    # evaluate (approximated) -H
    hhh <- lapply(1:length(Hmats), function(i){
      evalMat_row(Hmats[[i]], exp_dataList[[i]], par_val = par_val) %>% as.data.frame() %>%
        cbind(exp_dataList[[i]][subject_id], .)
    }) %>% Reduce("rbind", .) %>%
      aggregate( as.formula(paste(". ~", subject_id)), ., sum)
    
    log_detH <- lapply(1:n, function(i){
      h_i <- -(matrix(as.numeric(hhh[i, -1]), q, q) + negH3[[i]])
      log(det(h_i))
    }) %>% unlist() %>% sum()
    
    # evaluate h-likelihood
    lhlike <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]), invSIGMA, distribution, degree)
    # evaluate negative profile h-likelihood
    -(lhlike - 0.5 * log_detH)
    
  }
  
  ff <- function(xx, ...){
    # all parameter values
    par_val <- Vassign(names(param), xx) %>% c(., other_param)
    
    if (naive == FALSE){
      new_dataList[[3]]$fitted_CD4 <- eval_fn_row(Object_CD4$reg_equation, new_dataList[[3]], par_val, get_gradient = F)
      max_CD4 <- new_dataList[[3]] %>%
        group_by(!!sym(subject_id)) %>%
        filter(fitted_CD4 == max(fitted_CD4)) %>%
        dplyr::select(!!sym(subject_id), fitted_CD4) 
      colnames(max_CD4)[2] <- "maxCD4"
      new_dataList[[4]] <- new_dataList[[4]] %>% dplyr::select(-tidyselect::one_of(c("maxCD4")))
      new_dataList[[4]] <- merge(new_dataList[[4]], max_CD4) 
    } 
    
    # evaluate h-likelihood
    lhlike <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]), invSIGMA, distribution, degree)
    # evaluate negative profile h-likelihood
    - lhlike
    
  }
  
  
  ## gr() function is not used because it's slower than calculating numerical gradients
  if (adjust == TRUE){
    res <- try(optim_iterate_fix(str_val = unlist(param), ff_adj, gr = NULL, lower = lower, upper = upper,
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
                                 # dH_theta = dH_theta,
                                 # dfs = dfs,
                                 uniqueID = uniqueID, # fix_pars = fix_pars,
                                 all_glmeObjects = all_glmeObjects))
  } else{
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
                             # dH_theta = dH_theta,
                             # dfs = dfs,
                             uniqueID = uniqueID, # fix_pars = fix_pars,
                             all_glmeObjects = all_glmeObjects) 
  }
  res$gamma
}


