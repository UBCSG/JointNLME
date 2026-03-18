#' - d^2 l_h / db db^T
eval_H_b <- function(RespLog, random_effects, new_dataList,
                     Bi_df, invSIGMA, par_val,
                     distribution, degree, extra_pars = NULL){
  n <- nrow(Bi_df)
  q <- length(random_effects)

  # Derive H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars = c(random_effects, extra_pars))
  deno <- diag(as.matrix(Bi_df[, -1])%*%invSIGMA%*%t(Bi_df[, -1]))
  numo <- as.matrix(Bi_df[, -1])%*%invSIGMA
  # evaluate the -H matrix for the random effects part
  if(distribution == "normal"){
    negH3 <- -invSIGMA*n
  } else if (distribution == "t-dist"){
    negH3 <- lapply(1:n, function(i){
      numo1 <- (degree + q)*invSIGMA*(degree + deno[i]) -
        (degree + q)*as.matrix(numo[i, ])%*%(2*numo[i, ])
      -numo1/((degree + deno[i])^2)
    }) %>% Reduce('+', .)
  }
  negH <- lapply(1:length(Hmats), function(i){
    evalMat(Hmats[[i]], new_dataList[[i]], par_val = par_val)
  }) %>% Reduce("+", .)

  if(is.null(extra_pars)){
    -(negH + negH3)
  } else {
    -(negH + bdiag(negH3, diag(0, length(extra_pars))))
  }
}

#' d H_b / d theta
eval_dHb_dtheta <- function(RespLog, random_effects,
                            fixed_param, new_dataList, par_val,
                            extra_pars = NULL){
  q <- length(random_effects) + length(extra_pars)
  p <- length(fixed_param)

  Hmats <- getHmat(RespLog, pars = c(random_effects, extra_pars))
  dH_i <- array(0, dim = c(q, q, p))
  for(i in 1:p){
    dH_i[, , i] <- lapply(1:length(Hmats), function(k){
      lapply(Hmats[[k]], function(x){
        D(x, names(fixed_param)[i])
      }) %>%
        evalMat(new_dataList[[k]], par_val = par_val)
    }) %>% Reduce("+", .)
  }
  -dH_i
}

#' d^2 H_b / d theta d theta^T
eval_d2Hb_d2theta <- function(RespLog, random_effects,
                              fixed_param, new_dataList, par_val,
                              extra_pars = NULL){
  q <- length(random_effects) + length(extra_pars)
  p <- length(fixed_param)
  Hmats <- getHmat(RespLog, pars = c(random_effects, extra_pars))
  dH_ij <- array(0, dim = c(q, q, p, p))
  for(i in 1:p){
    for(j in i:p){
      Dij <- lapply(1:length(Hmats), function(k){
        lapply(Hmats[[k]], function(x){
          D(x, names(fixed_param)[i]) %>%
            D(names(fixed_param)[j]) }) %>%
          evalMat(new_dataList[[k]], par_val = par_val)
      }) %>% Reduce("+", .)
      dH_ij[,,i,j] <- dH_ij[,,j,i] <- Dij
    }
  }
  -dH_ij
}

#' Estimate cov-matrix of fixed parameters based on profile h-likelihood
#'
#' cov(theta) = d^2 l_p /d theta d theta^T
get_sd_phlike <- function(fixed_param, random_effects,
                          other_param,
                          RespLog, dataList,
                          Bi_df, invSIGMA,
                          # lambda = 0,
                          distribution = "normal", degree = 0,
                          extra_pars = NULL){

  subject_id <- names(Bi_df)[1]
  random_effects <- names(Bi_df)[-1]
  q <- length(random_effects)
  n <- nrow(Bi_df)
  p <- length(fixed_param)
  new_dataList <- lapply(dataList, function(dat){dplyr::left_join(dat, Bi_df, by = subject_id)})
  par_val <- c(fixed_param, other_param)

  # part1: d^2 l_h / d theta d theta^T
  Hmats_theta <- getHmat(RespLog, pars = names(fixed_param))
  Hval_theta <- lapply(1:length(Hmats_theta), function(i){
    evalMat(Hmats_theta[[i]], new_dataList[[i]], par_val = par_val)
  }) %>% Reduce("+", .)

  # part2: d^2 log(det(H_b)) =
  #  (d^2 |H_b|/d theta d theta^T) / |H_b| -
  # (d |H_b|/d theta)*(d |H_b|/d theta^T)/(|H_b|^2)

  # extra_pars <- names(disp_param)[5:9]
  # extra_pars <- NULL
  Hval_b <- eval_H_b(RespLog, random_effects, new_dataList, Bi_df,
                     invSIGMA, par_val, distribution, degree, extra_pars)
  invH <- solve(Hval_b)
  dH_i <- eval_dHb_dtheta(RespLog, random_effects, fixed_param, new_dataList, par_val, extra_pars)
  d2H_ij <- eval_d2Hb_d2theta(RespLog, random_effects, fixed_param, new_dataList, par_val, extra_pars)

  Hessian2 <- matrix(NA, p, p)
  for(i in 1:p){
    for(j in i:p){
      mat_i <- invH%*%dH_i[,,i]
      mat_j <- invH%*%dH_i[,,j]
      Hessian2[i, j] <- Hessian2[j, i] <-
        matrix.trace(as.matrix(- mat_i%*%mat_j + invH%*%d2H_ij[,,i,j]))
    }
  }

  Hessian <- Hval_theta - 0.5*Hessian2
  covMat <- solve(-Hessian)
  sqrt(diag(covMat))
}
