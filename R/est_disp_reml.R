#' #' Estimate the dispersion parameters in the joint models,
#' #' by maximizing the adjusted profile h-likelihood function.
#' #'
#' #' @importFrom Matrix bdiag
#' #' @importFrom dplyr select
#' #' @importFrom matrixcalc matrix.trace
#' #' @importFrom dplyr left_join
est_disp_reml <- function(param, other_param,
                         RespLog, # deriv_vars,
                         new_dataList, # dataList,
                         Bi_df, invSIGMA, last_decayobs, T_end,
                         lower = -Inf, upper = Inf,
                         distribution = 'normal', degree = 0, correlation = NULL,
                         cov_method = c("spherical", "cholesky"),
                         # disp_pars, fix_pars, 
                         Hmats = Hmats,
                         Hmats3 = Hmats3,
                         REML, all_glmeObjects,
                         Silent = T){

  # q <- length(deriv_vars)
  q1 <- length(param)
  q2 <- as.integer(ncol(Bi_df) - 1)
  q3 <- length(other_param)
  n <- nrow(Bi_df)
  subject_id <- names(Bi_df)[1]
  uniqueID <- unique(Bi_df[, 1])
  
  # if(REML == FALSE){
  #   Hmats <- fix_pars$Hmats
  # } else {
  #   Hmats <- disp_pars$Hmats
  #   dhlike <- disp_pars$dhlike
  #   dHmats <- disp_pars$dHmats
  # }

  # Return matrices of strings for cov(bi)
  if(q2 > 1 | cov_method == "cholesky"){
    L2 <- create_covb_param(q2, cov_method, correlation)
    q0 <- length(L2$Mpar)
  } else {
    stop("When there is only 1 random effect in the model, must use cholesky method for estimating var-cov matrix of random effects.
         Change pre_memList format accordingly.")
  }

  # Return the negative value of the adjusted profile h-likelihood
  ff <- function(xx, ...){
    if(!is.null(param)){
      # The first q1 are dispersion parameters in the models,
      # and the rests are for the cov.matrix of random effects
      par_val <- c(Vassign(names(param), xx[1:q1]), other_param)
      Lval <- Vassign(L2$Mpar, xx[-(1:q1)])
    } else {
      # No other dispersion parameters except for the covariance
      # matrix of random effects
      par_val <- other_param
      Lval <- Vassign(L2$Mpar, xx)
    }

    if(distribution == 'normal'){
      if (q0 == 0){
        mat <- invSIGMA
      } else{
        invmat <- evalMat(as.list(L2$M), par_val = Lval) # i.e. SIGMA
        mat <- solve(invmat)
      }
    } else if(distribution == 't-dist'){
      invmat <- evalMat(as.list(L2$M), par_val = Lval)*(degree - 2)/degree
      mat <- solve(invmat)
    }

    # component by random effect
    if(distribution == 'normal'){
      nD3 <- lapply(1:n, function(k){
         # - mat
        Matrix::bdiag(-mat, diag(0, q3 + 1, q3 + 1))
        })
    } else if(distribution == 't-dist'){
      numo <- as.matrix(Bi_df[, -1]) %*% mat
      deno <- numo%*%t(Bi_df[, -1]) %>% diag()
      nD3 <- list()
      for(i in 1:n){
        numo1 <- (degree + q2)*mat*(degree + deno[i]) - (degree + q2)*as.matrix(numo[i, ])%*%(2*numo[i, ])
        nD3_i <- -numo1/((degree + deno[i])^2)
        nD3[[i]] <- nD3_i
      }
    }
    # if(REML == TRUE){
    #   nD3 <- lapply(nD3, function(X){ Matrix::bdiag( X, diag(0, q-q2, q-q2)) %>% as.matrix() })
    # }


    # To approximate -H by E(-H|b) for NLMEs
    # Ref: Maengseok Noh and Youngjo Lee, 2008, Hierarchical-likelihood approach for nonlinear mixed-effects models
    # exp_dataList <- dataList
    # exp_dataList <- lapply(dataList, function(dat){
    #   dplyr::left_join(dat, as.data.frame(Bi_df), by = subject_id)
    # })
    # 
    # exp_dataList[[4]][all_glmeObjects[[4]]$response] <- exp(eval_fn_row(all_glmeObjects[[4]]$reg_equation, exp_dataList[[4]], par_val, get_gradient = F))
    # Ti_new <- exp_dataList[[4]] %>% dplyr::select(ID, time)
    # names(Ti_new)[2] <- "Ti_new"
    # exp_dataList[1:3] <- lapply(exp_dataList[1:3], function(dat){
    #   dat <- dplyr::left_join(dat, Ti_new, by = subject_id)
    # })
    # exp_dataList[[1]] <- exp_dataList[[1]] %>%
    #   filter(time_decay <= Ti_new + last_decayobs) %>%
    #   mutate(y_decay_old = y_decay) %>%
    #   mutate(censor_decay = 0)
    # exp_dataList[[2]] <- exp_dataList[[2]] %>%
    #   filter(time_rebound > Ti_new + last_decayobs) %>%
    #   mutate(time_rebound = time_rebound - (Ti_new + last_decayobs)) %>%
    #   mutate(y_rebound_old = y_rebound) %>%
    #   mutate(censor_rebound = 0)
    # exp_dataList[[3]] <- exp_dataList[[3]] %>%
    #   filter(time_CD4 <= Ti_new + last_decayobs) %>%
    #   mutate(CD4_old = CD4)
    # 
    # for(j in 1:(length(all_glmeObjects) - 1)){
    #   exp_dataList[[j]][all_glmeObjects[[j]]$response] <-
    #     eval_fn_row(all_glmeObjects[[j]]$reg_equation, exp_dataList[[j]], par_val, get_gradient = F)
    # }
    # 
    # evaluate the H matrix
    hhh <- lapply(1:length(Hmats3), function(k){
      evalMat_row(Hmats3[[k]], new_dataList[[k]], par_val = par_val) %>% as.data.frame() %>%
        cbind(new_dataList[[k]][subject_id], .)
    }) %>% Reduce("rbind", .) %>%
      aggregate( as.formula(paste(". ~", subject_id)), ., sum)

    # if(REML == FALSE){
    #   log_det <- lapply(1:n, function(k){
    #     h_i <- -(matrix(as.numeric(hhh[k, -1]), q2 + 1, q2 + 1) + nD3[[k]])/(2*pi)
    #     log(det(nearPD(h_i)$mat))
    #   }) %>% unlist() %>% sum()
    # } else {
    #   log_det <- get_log_det_disp_cpp(as.matrix(hhh[, -1]), nD3, q2, q-q2)
    # }

    log_det <- lapply(1:n, function(k){
      h_i <- -(matrix(as.numeric(hhh[k, -1]), q2 + 1 + q3, q2 + 1 + q3) + nD3[[k]])/(2*pi)
      log(det(nearPD(h_i)$mat))
    }) %>% unlist() %>% sum()
    
    lhlike <- numeric(1)

    # evaluate the h-likelihood
    lhlike <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]),
                             mat, distribution, degree)
    # evaluate the negative adjusted profile h-likelihood
    -(lhlike - 0.5*log_det) # - 0.5*(n*q2 + q-q2)*log(2*pi)

    # lhlike <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]), mat, distribution, degree)
    # - lhlike
  }


  if(cov_method == "cholesky"){
    lower_new <- c(lower, rep(-Inf, q0))
    upper_new <- c(upper, rep(Inf, q0))
  } else if(cov_method == "spherical"){
    lower_new <- c(lower, rep(0, q0))
    upper_new <- c(upper, rep(pi, q0))
  }

  Lval0 <- initiate_covb_par(L2, invSIGMA, cov_method, correlation)
  str_val <- c(param, Lval0) %>% unlist()
  
  # res <- optim_iterate_fix(str_val, ff, gr = NULL, lower = lower_new, upper = upper_new,
  #                             check = 1 - as.numeric(Silent), 
  #                             # extra info
  #                             param = param, other_param = other_param,
  #                             Hmats = Hmats, new_dataList = new_dataList,
  #                          # dataList = datalist,
  #                             RespLog = RespLog, Bi_df = Bi_df,
  #                             distribution = distribution, degree = degree,
  #                             # dhlike = dhlike,
  #                             # dHmats = dHmats,
  #                             uniqueID = uniqueID,
  #                             q = q, q1 = q1, q2 = q2, q0 = q0, n = n, L2 = L2,
  #                             all_glmeObjects = all_glmeObjects)
                           #disp_pars = disp_pars, fix_pars = fix_pars)
  res <- optim_iterate_nlminb(str_val, ff, gr = NULL, lower = lower_new, upper = upper_new,
                              check = 1 - as.numeric(Silent), #gr_tol = 1,
                              rel_tol = 1e-3, x_tol = 1e-3,
                              # rel_tol = 1e-1, x_tol = 1e-1,
                              step_min = 1e-5, step_max = 1,
                              # abs.tol = 1e-1,
                              # extra info
                              param = param, other_param = other_param,
                              Hmats = Hmats, Hmats3 = Hmats3, new_dataList = new_dataList,
                              RespLog = RespLog, Bi_df = Bi_df,
                              distribution = distribution, degree = degree,
                              # dhlike = dhlike,
                              # dHmats = dHmats,
                              # uniqueID = uniqueID,
                              q = q, q1 = q1, q2 = q2, q0 = q0, n = n, L2 = L2,
                              all_glmeObjects = all_glmeObjects)
  

  disp_param <- res$gamma[1:q1]

  if(distribution == 'normal'){
    if (q0 == 0){
      invmat = solve(invSIGMA)
    } else{
      invmat <- evalMat(as.list(L2$M), par_val = res$gamma[-(1:q1)]) 
    }
  } else if (distribution == 't-dist'){
    invmat <- evalMat(as.list(L2$M), par_val = res$gamma[-(1:q1)])*(degree - 2)/degree
  }
  mat <- solve(invmat)

  list(disp_param = disp_param, invSIGMA = mat, 
       Lval = res$gamma[-(1:q1)])
}

