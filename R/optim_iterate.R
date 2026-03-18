#' General-purpose optimization
#'
#' @param str_val a named vector of initial values
#' @param ff a function to be minimized
#' @param gr gradient function of fn
#' @param lower,upper bounds on the parameters to be estimated
#' @param check 0 or 1; If 1, print information on the progress of the optimization
#' @param gr_tol tolerance level of final gradient
#' @importFrom optimx optimx
#' @importFrom magrittr %>%
optim_iterate_rand <- function(str_val, ff, par_val, gr = NULL, lower = -Inf, upper = Inf, check = 1,
                          gr_tol = 1, random_start = TRUE, methods = NULL, ...){

  if(is.null(methods)){
    if(all(lower == -Inf) & all(upper == Inf)){
      methods <- c("nlminb", "BFGS")
    } else {
      methods <- c("nlminb", "L-BFGS-B")
    }
  }

  converged <- FALSE
  M <- 0
  str_val0 <- str_val


  while(!converged & M < 100){
    error_mess <- "try-error"
    for(method in methods){
      if("try-error" %in% error_mess | !converged){
        # cat(method, "is used.", '\n')
        result <- try(optimx::optimx(par = str_val0, fn = ff, gr = gr,
                                     lower = lower, upper = upper, method = method,
                                     itnmax = 100,
                                     control = list(trace = 1 - check, kkt = FALSE,
                                                    follow.on = FALSE, starttests = FALSE)), silent = T)
        error_mess <- attr(result, "class")

        if("try-error" %in% error_mess){
          print(result)
          converged <- FALSE
        } else {
          if(result$convcode != 0){
            converged <- FALSE
          } else {
            str_val2 <- as.numeric(result[1:length(str_val)])
            converged=TRUE
          }
        }
      }
    }

    str_bi <- str_val[-length(str_val)]
    lower_bi <- lower[-length(lower)]
    upper_bi <- upper[-length(upper)]
    str_Ti <- str_val[length(str_val)]
    lower_Ti <- lower[length(lower)]
    upper_Ti <- upper[length(upper)]
    if(converged == FALSE){
      if(random_start == TRUE){
        str_bi0 <- rnorm(length(str_bi), 0, 1) %>%
          pmax(lower_bi) %>% pmin(upper_bi)
        str_Ti0 = rweibull(1, shape = par_val$Wshape, scale = par_val$Wscale * exp(par_val$gamma1 * str_bi0[2] + par_val$gamma2 * str_bi0[6]))%>%
          pmax(lower_Ti) %>% pmin(upper_Ti)
        str_val0 <- c(str_bi0, str_Ti0)
      } else {
        str_bi0 <- str_bi + runif(length(str_bi), - abs(str_bi/10), abs(str_bi/10))
        str_Ti0 <- str_Ti + runif(1, - str_Ti/10,  str_Ti/10)
        str_bi0 <- str_bi0 %>% pmax(lower_bi) %>% pmin(upper_bi)
        str_Ti0 <- str_Ti0 %>% pmax(lower_Ti) %>% pmin(upper_Ti)
      }
    }
    M <- M + 1
  }

  if(converged){
    names(str_val2) <- names(str_val)
    list(gamma = str_val2, fval = result$value)
  } else {
    stop("Failed to converge ...")
    # list(gamma = rep(NA, length(str_val)), fval = NA)
  }
}

optim_iterate_fix <- function(str_val, ff, gr = NULL, lower = -Inf, upper = Inf, check = 1,
                          gr_tol = 1, random_start = FALSE, methods = NULL, ...){
  
  if(is.null(methods)){
    if(all(lower == -Inf) & all(upper == Inf)){
      methods <- c("nlminb", "BFGS")
    } else {
      methods <- c("nlminb", "L-BFGS-B")
    }
  }
  
  converged <- FALSE
  M <- 0
  str_val0 <- str_val
  # print(ff(str_val0))
  # print(gr(str_val0))
  
  while( !converged & M < 100){
    error_mess <- "try-error"
    for(method in methods){
      if("try-error" %in% error_mess | !converged){
        # cat(method, "is used.", '\n')
        result <- try(optimx::optimx(par = str_val0, fn = ff, gr = gr,
                                     lower = lower, upper = upper, method = method,
                                     itnmax = 500,
                                     control = list(trace = check, kkt = FALSE,
                                                    follow.on = FALSE, starttests = FALSE)), silent = T)
        
         error_mess <- attr(result, "class")
        
        if("try-error" %in% error_mess){
          print(result)
          converged <- FALSE
        } else {
          if(result$convcode != 0){
            converged <- FALSE
          } else {
            str_val2 <- as.numeric(result[1:length(str_val)])
            # if(is.null(gr)){
            #   final_gr <- num_gradient(str_val2, ff, ...)
            #   print(final_gr)
            # } else{
            #   final_gr <- gr(str_val2, ...)
            # }
            # converged <- ifelse(sum(is.na(final_gr)) > 0, FALSE, max(abs(final_gr)) < gr_tol)
            converged=TRUE
          }
        }
      }
    }
    
    if(converged == FALSE){
      if(random_start == TRUE){
        str_val0 <- rnorm(length(str_val), 0, 1) %>%
          pmax(lower) %>% pmin(upper)
      } else {
        str_val0 <- str_val + runif(length(str_val), - abs(str_val/10), abs(str_val/10))
        str_val0 <- str_val0 %>% pmax(lower) %>% pmin(upper)
      }
    }
    M <- M +1
  }
  
  if(converged){
    names(str_val2) <- names(str_val)
    list(gamma = str_val2, fval = result$value)
  } else {
    stop("Failed to converge ...")
  }
}


#' Optimization using nlminb
#'
#' @param str_val a named vector of initial values
#' @param ff a function to be minimized
#' @param gr gradient function of fn
#' @param lower,upper bounds on the parameters to be estimated
#' @param check 0 or 1; If 1, print information on the progress of the optimization
#' @param gr_tol tolerance level of final gradient
#' @importFrom optimx optimx
#' @importFrom magrittr %>%
optim_iterate_nlminb <- function(str_val, ff, gr = NULL, lower = -Inf, upper = Inf, check = 1,
                                 gr_tol = 1, random_start = FALSE,
                                 rel_tol = 1e-10, x_tol = 1e-10,# xf_tol = 1e-3,
                                 step_min = 1e-5, step_max = 1, iter_max = 500,
                                 ...){
  converged <- FALSE
  M <- 0
  str_val0 <- str_val

  while( !converged & M < 10){
    result <- try(nlminb(start = str_val0, objective = ff, gradient = gr, 
                         control = list(rel.tol = rel_tol, x.tol = x_tol,
                                        # abs.tol = 1e-1,
                                        # xf.tol = xf_tol,
                                        step.min = step_min,
                                        step.max = step_max,
                                        iter.max = iter_max,
                                        trace = check),
                         lower = lower, upper = upper), silent = T)
    error_mess <- attr(result, "class")
    # print(result)


    if("try-error" %in% error_mess){
      converged <- FALSE
      print(result)
    } else {
      # print(result)
      # print(num_gradient(result$par, ff, ...))
      # print(final_gr <- gr(result$par, ...))
      if(result$convergence != 0){
        converged <- FALSE
      } else {
        str_val2 <- result$par

        # if(is.null(gr)){
        #   final_gr <- num_gradient(str_val2, ff, ...)
        #   print(final_gr)
        # } else{
        #   final_gr <- gr(str_val2, ...)
        # }
        # converged <- ifelse(sum(is.na(final_gr)) > 0, FALSE, max(abs(final_gr)) < gr_tol)
        converged <- TRUE
      }
    }
    if(converged == FALSE){
      if(random_start == TRUE){
        str_val0 <- rnorm(length(str_val), 0, 1) %>%
          pmax(lower) %>% pmin(upper)
      } else {
        str_val0 <- str_val + runif(length(str_val), - abs(str_val/10), abs(str_val/10))
        str_val0 <- as.numeric(str_val0) %>% pmax(lower) %>% pmin(upper)
      }
    }
    M <- M +1
  }

  if(converged){
    names(str_val2) <- names(str_val)
    list(gamma = str_val2, fval = result$value)#, gval = final_gr)
  } else {
    stop("Failed to converge ...")
  }
}

#' Optimization using L-BFGS
#'
#' @param str_val a named vector of initial values
#' @param ff a function to be minimized
#' @param gr gradient function of fn
#' @param lower,upper bounds on the parameters to be estimated
#' @param check 0 or 1; If 1, print information on the progress of the optimization
#' @param gr_tol tolerance level of final gradient
#' @importFrom optimx optimx
#' @importFrom magrittr %>%
optim_iterate_lbfgs <- function(str_val, ff, gr, lower = -Inf, upper = Inf, check = 1,
                                 gr_tol = 0.1, random_start = FALSE, ...){
  converged <- FALSE
  M <- 0
  str_val0 <- str_val

  while( !converged & M < 50){
    result <- try(lbfgsb3c::lbfgsb3c(par = str_val0, fn = ff, gr = gr,
                                     lower = lower, upper = upper,
                                     control = list(trace = 1 - check),
                                     ...), silent = T)

    # result <- try(optim(par = str_val0, fn = ff, gr = gr,
    #                     method = "L-BFGS-B",
    #                     lower = lower, upper = upper, hessian = FALSE,
    #                     control = list(trace = 1, REPORT = 1)), silent = T)

    error_mess <- attr(result, "class")
    # print(error_mess)
    # print(result)
    # print(final_gr <- gr(result$par, ...))

    if("try-error" %in% error_mess){
      converged <- FALSE
    } else {
      if(result$convergence != 0){
        converged <- FALSE
      } else {
        str_val2 <- result$par
        # if(is.null(gr)){
        #   final_gr <- num_gradient(str_val2, ff, ...)
        #   print(final_gr)
        # } else{
        #   final_gr <- gr(str_val2, ...)
        # }
        # converged <- ifelse(sum(is.na(final_gr)) > 0, FALSE, max(abs(final_gr)) < gr_tol)
        converged <- TRUE
      }
    }
    if(converged == FALSE){
      if(random_start == TRUE){
        str_val0 <- rnorm(length(str_val), 0, 1) %>%
          pmax(lower) %>% pmin(upper)
      } else {
        str_val0 <- str_val + runif(length(str_val), - abs(str_val/10), abs(str_val/10))
        str_val0 <- str_val0 %>% pmax(lower) %>% pmin(upper)
      }
    }
    M <- M +1
  }

  if(converged){
    names(str_val2) <- names(str_val)
    list(gamma = str_val2, fval = result$value, gval = final_gr)
  } else {
    stop("Failed to converge ...")
  }
}

# TODO: May consider parallel version
# Note: Using optimParallel with n parallel processes increases the
# memory usage by about factor n compared to a call to optim. If the
# memory limit is reached this may severely slowdown the optimization.
# cl <- parallel::makeCluster(3)
# parallel::setDefaultCluster(cl = cl)
# parallel::clusterEvalQ(cl, devtools::load_all())
# parallel::clusterExport(cl, c("RespLog", "sub_long", "sub_surv",
#                               "random_effects", "invSIGMA", "par_val",
#                               "distribution", "degree", "q",
#                               "gr_long", "gr_surv"))
# result <- try(optimParallel::optimParallel(par = bval, fn = ff, gr = gr,
#                     control = list(maxit = 100, trace = 0),
#                     parallel=list(forward = TRUE, loginfo = TRUE)), silent = T)
# parallel::setDefaultCluster(cl = NULL)
# parallel::stopCluster(cl)
