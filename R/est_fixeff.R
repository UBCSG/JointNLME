#' Estimate fixed parameters by maximizing the profile h-likelihood function
#'
#' @importFrom dplyr select left_join
#' @importFrom Deriv Simplify
#' @importFrom matrixcalc matrix.trace
est_fixeff <- function(param, other_param,
                       RespLog, dataList,
                       Bi_df, invSIGMA,
                       lower = -Inf, upper = Inf,
                       distribution = "normal", degree = 0,
                       fix_pars,
                       Silent = TRUE){

  subject_id <- names(Bi_df)[1]
  random_effects <- names(Bi_df)[-1]
  q <- length(random_effects)
  n <- nrow(Bi_df)
  p <- length(param)
  new_dataList <- lapply(dataList, function(dat){dplyr::left_join(dat, Bi_df, by = subject_id)})
  Hmats <- fix_pars$Hmats
  dH_theta <- fix_pars$dH_theta
  dfs <- fix_pars$dfs

  if(distribution == "normal"){
    negH3 <- -invSIGMA*n
  } else if (distribution == "t-dist"){
    numo <- as.matrix(Bi_df[, -1]) %*% invSIGMA
    deno <- numo %*% t(Bi_df[, -1]) %>% diag()
    negH3 <- lapply(1:n, function(i){
      numo1 <- (degree + q)*invSIGMA*(degree + deno[i]) -
        (degree + q)*as.matrix(numo[i, ])%*%(2*numo[i, ])
      -numo1/((degree + deno[i])^2)
    }) %>% Reduce('+', .)
  }

  # Negative adjusted profile h-likelihood with param as input
  ff <- function(xx, ...){
    # all parameter values
    par_val <- Vassign(names(param), xx) %>% c(., other_param)

    # evaluate the -H matrix
    negH <- lapply(1:length(Hmats), function(i){
      evalMat(Hmats[[i]], new_dataList[[i]], par_val)
    }) %>% Reduce("+", .)
    H <- - (negH + negH3)

    # evaluate h-likelihood
    lhlike <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]), invSIGMA, distribution, degree)
    # evaluate negative profile h-likelihood
    -(lhlike - 0.5*log(det(H/2/pi)))
  }

  # Gradient of ff()
  gr <- function(xx, ...){
    # all parameter values
    par_val <- Vassign(names(param), xx) %>% c(., other_param)

    # evaluate inverse -H matrix
    negH <- lapply(1:length(Hmats), function(i){
      evalMat(Hmats[[i]], new_dataList[[i]], par_val)
    }) %>% Reduce("+", .)
    H <- - (negH + negH3)
    invH <- solve(H)

    # evaluate dH/d theta
    fy <- lapply(1:length(param), function(i){
      dH_val <- matrix(0, q, q)
      df_val <- 0
      for(j in 1:length(RespLog)){
        dH_val_j <- evalMat(dH_theta[[j]][[i]], new_dataList[[j]], par_val)
        dH_val <- dH_val - dH_val_j # Note: Hmats returns -H matrix
        df_val_j <- eval_fn(dfs[[j]][[i]], new_dataList[[j]], par_val, get_gradient = FALSE)
        df_val <- df_val + sum(df_val_j)
      }
      -(df_val - 0.5*matrix.trace(invH %*% dH_val))
    }) %>% unlist()

    fy
  }

  # check_gradient(unlist(param), ff, gr) %>% print
  res <- optim_iterate(unlist(param), ff, gr, lower = lower, upper = upper,
                       check = as.numeric(Silent), gr_tol = 1,
                       # extra inputs
                       param = param, other_param = other_param,
                       RespLog = RespLog,
                       Bi_df = Bi_df, invSIGMA = invSIGMA,
                       distribution = distribution, degree = degree,
                       random_effects = random_effects,
                       new_dataList = new_dataList,
                       Hmats = Hmats, negH3 = negH3,
                       dH_theta = dH_theta,
                       dfs = dfs)
  res
}
