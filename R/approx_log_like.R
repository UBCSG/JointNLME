#' Evaluate the approximate log marginal likelihood function
#'
approx_log_like <- function(RespLog, par_val, new_dataList, Bi_df,invSIGMA, 
                            distribution, degree,
                            all_glmeObjects, fix_pars, naive){
  par_val <- as.list(par_val)
  random_effects <- names(Bi_df)[-1]
  q <- length(random_effects)
  n <- nrow(Bi_df)
  subject_id <- names(Bi_df)[1]
  uniqueID <- Bi_df[, 1]
  Hmats <- fix_pars$Hmats
  
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
        eval_fn_row(all_glmeObjects[[j]]$reg_equation, new_dataList[[j]], par_val, get_gradient = F)
    }
    exp_dataList[[4]][all_glmeObjects[[4]]$response] <- exp(exp_dataList[[4]][all_glmeObjects[[4]]$response])
  } else{
    exp_dataList <- new_dataList
    for(j in 1:length(all_glmeObjects)){
      exp_dataList[[j]][all_glmeObjects[[j]]$response] <-
        eval_fn_row(all_glmeObjects[[j]]$reg_equation, new_dataList[[j]], par_val, get_gradient = F)
    }
    exp_dataList[[3]][all_glmeObjects[[3]]$response] <- exp(exp_dataList[[3]][all_glmeObjects[[3]]$response])
  }
  
  hhh <- lapply(1:length(Hmats), function(i){
    evalMat_row(Hmats[[i]], exp_dataList[[i]], par_val) %>% as.data.frame() %>%
      cbind(exp_dataList[[i]][subject_id], .)
  }) %>% Reduce("rbind", .) %>%
    aggregate( as.formula(paste(". ~", subject_id)), ., sum)
  
  Hval <- lapply(1:n, function(k){
    -(matrix(as.numeric(hhh[k, -1]), q, q) - invSIGMA)
  }) %>% Reduce(bdiag, .) %>% as.matrix()
  
  hloglike_value <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]), invSIGMA, distribution, degree)
  
  hloglike_value - 0.5*log(det(Hval/2/pi))
}

#' Evaluate log h-likelihood
#'
eval_log_hlike <- function(RespLog, new_sub_dataList, par_val, Bi,
                           invSIGMA, distribution, degree){
  stopifnot(class(Bi)[1] == "matrix")
  if(class(par_val) != "list"){
    par_val <- as.list(par_val)
  }
  q <- ncol(Bi)

  lhlike <- sapply(1:length(RespLog), function(i){
    if (nrow(new_sub_dataList[[i]]) != 0){
      eval_fn(RespLog[[i]], new_sub_dataList[[i]], par_val, get_gradient = FALSE)
    }
  })
  deno <- diag(Bi%*%invSIGMA%*%t(Bi))
  stopifnot(length(deno) == nrow(Bi))

  lhlike3 <- 0.5*log(det(invSIGMA)) - 0.5*deno - q/2*log(2*pi)

  sum(unlist(lhlike)) + sum(lhlike3)
}


