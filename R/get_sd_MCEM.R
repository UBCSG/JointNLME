#' Estimate s.e. of non-robust fixed parameters,
#' using MCEM approximation
#' 
get_sd_MCEM <- function(fixed_param, RespLog, pre_data, post_data,
                        Bis, random_effects, subject_id, disp_param,
                        distribution, degree, Silent){
  
  q <- length(random_effects)
  p <- length(fixed_param)
  n <- nrow(post_data)
  gr.long <- deriv(formula(paste("~", RespLog[[1]])), names(fixed_param), hessian = T)
  gr.surv <- deriv(formula(paste("~", RespLog[[2]])), names(fixed_param), hessian = T)
  
  gr <- function(xx, subB_i){
    fy <- numeric(p)
    # assign values to parameters
    par_val <- c(Vassign(names(fixed_param), xx), disp_param)
    subB <- dplyr::select(pre_data, subject_id) %>%
      dplyr::left_join(subB_i, by = subject_id)
    post_subB <- dplyr::select(post_data, subject_id) %>%
      dplyr::left_join(subB_i, by = subject_id)
    fy1 <- with(pre_data, with(par_val, with(subB, attr(eval(gr.long), "gradient"))))
    fy2 <- with(post_data, with(par_val, with(post_subB, attr(eval(gr.surv), "gradient"))))
    score_i <- apply(fy1, 2, function(x){
      tapply(x, pre_data[, subject_id], sum)
    }) + apply(fy2, 2, function(x){
      tapply(x, post_data[, subject_id], sum)
    })
    
    jac1 <- with(pre_data, with(par_val, with(subB, attr(eval(gr.long), "hessian"))))
    jac2 <- with(post_data, with(par_val, with(post_subB, attr(eval(gr.surv), "hessian"))))
    jac_i <- apply(jac1, MARGIN = c(2, 3), sum) + 
      apply(jac2, MARGIN = c(2, 3), sum)
    
    # score_i%*%score_i
    return(list(score = score_i, jac = jac_i))
  }
  
  mcmc_size <- length(Bis)
  mat1 <- mat2 <- matrix(0, p, p)
  score <- rep(0, p)
  for(i in 1:mcmc_size){
    out <- gr(unlist(fixed_param), Bis[[i]])
    mat1 <- mat1 + t(out$score)%*%out$score
    mat2 <- mat2 + out$jac
    score <- score + out$score
  }
  
  result <- (mat1 + mat2 - t(score)%*%score/mcmc_size)/mcmc_size
  sqrt(diag(-solve(result)))
}