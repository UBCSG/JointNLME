#' Estimate the dispersion parameters in the joint models,
#' by maximizing the adjusted profile h-likelihood function.
#'
#' @importFrom Matrix bdiag
#' @importFrom dplyr select
#' @importFrom matrixcalc matrix.trace
#' @importFrom dplyr left_join
est_disp <- function(param, other_param,
                     RespLog, deriv_vars,
                     dataList,
                     Bi_df, invSIGMA, lower = -Inf, upper = Inf,
                     distribution = 'normal', degree = 0, correlation = NULL,
                     cov_method = c("spherical", "cholesky"),
                     disp_pars = disp_pars,
                     Silent = T){

  q <- length(deriv_vars)
  q1 <- length(param)
  q2 <- ncol(Bi_df) - 1
  n <- nrow(Bi_df)
  subject_id <- names(Bi_df)[1]
  new_dataList <- lapply(dataList, function(dat){dplyr::left_join(dat, Bi_df, by = subject_id)})

  dhlike <- disp_pars$dhlike
  Hmats <- disp_pars$Hmats
  dHmats <- disp_pars$dHmats

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
      # and the rests are for the cov. matrix of random effects
      par_val <- c(Vassign(names(param), xx[1:q1]), other_param)
      Lval <- Vassign(L2$Mpar, xx[-(1:q1)])
    } else {
      # No other dispersion parameters except for the covariance
      # matrix of random effects
      par_val <- other_param
      Lval <- Vassign(L2$Mpar, xx)
    }
    if(distribution == 'normal'){
      invmat <- evalMat(as.list(L2$M), par_val = Lval) # i.e. SIGMA
      mat <- solve(invmat)
    } else if(distribution == 't-dist'){
      invmat <- evalMat(as.list(L2$M), par_val = Lval)*(degree - 2)/degree
      mat <- solve(invmat)
    }

    # evaluate the H matrix
    nD <- lapply(1:length(Hmats), function(i){
      evalMat(Hmats[[i]], new_dataList[[i]], par_val)
    }) %>% Reduce("+", .)
    # part 3: component by random effect
    if(distribution == 'normal'){
      nD3 <- Matrix::bdiag( -mat*n, diag(0, q-q2, q-q2))
    } else if(distribution == 't-dist'){
      numo <- as.matrix(Bi_df[, -1]) %*% mat
      deno <- numo%*%t(Bi_df[, -1]) %>% diag()
      nD3 <- matrix(0, nrow = q2, ncol = q2)
      for(i in 1:n){
        numo1 <- (degree + q2)*mat*(degree + deno[i]) - (degree + q2)*as.matrix(numo[i, ])%*%(2*numo[i, ])
        nD3_i <- -numo1/((degree + deno[i])^2)
        nD3 <- nD3 + nD3_i
      }
      nD3 <- Matrix::bdiag(nD3, diag(0, q-q2, q-q2))
    }
    H <- -as.matrix(nD + nD3)

    # evaluate the h-likelihood
    lhlike <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]),
                             mat, distribution, degree)
    # evaluate the negative adjusted profile h-likelihood
    -(lhlike - 0.5*log(det(H/2/pi)))
  }

  # Return the gradient of ff
  gr <- function(xx, ...){
    if(!is.null(param)){
      # The first q1 are dispersion parameters in the models,
      # and the rests are for the cov. matrix of random effects
      par_val <- c(Vassign(names(param), xx[1:q1]), other_param)
      Lval <- Vassign(L2$Mpar, xx[-(1:q1)])
    } else {
      # No other dispersion parameters except for the covariance
      # matrix of random effects
      par_val <- other_param
      Lval <- Vassign(L2$Mpar, xx)
    }
    if(distribution == 'normal'){
      invmat <- evalMat(as.list(L2$M), par_val = Lval) # i.e. SIGMA
      mat <- solve(invmat)
    } else if(distribution == 't-dist'){
      invmat <- evalMat(as.list(L2$M), par_val = Lval)*(degree - 2)/degree
      mat <- solve(invmat)
    }

    # evaluate the H matrix
    nD <- lapply(1:length(Hmats), function(i){
      evalMat(Hmats[[i]], new_dataList[[i]], par_val)
    }) %>% Reduce("+", .)

    if(distribution == 'normal'){
      nD3 <- Matrix::bdiag( -mat*n, diag(0, q-q2, q-q2))
    } else if(distribution == 't-dist'){
      numo <- as.matrix(Bi_df[, -1]) %*% mat
      deno <- numo %*% t(Bi_df[, -1]) %>% diag()

      nD3 <- matrix(0, nrow = q2, ncol = q2)
      for(i in 1:n){
        numo1 <- (degree + q2)*mat*(degree + deno[i]) - (degree + q2)*as.matrix(numo[i, ])%*%(2*numo[i, ])
        nD3_i <- -numo1/((degree + deno[i])^2)
        nD3 <- nD3 + nD3_i
      }
      nD3 <- Matrix::bdiag(nD3, diag(0, q-q2, q-q2))
    }
    H <- -as.matrix(nD + nD3)
    invH <- solve(H)

    # derivative w.r.t. dispersion parameters
    dh_val <- lapply(1:length(dhlike), function(i){
      eval_fn(dhlike[[i]], new_dataList[[i]], par_val, get_gradient = TRUE)
    }) %>% Reduce("+", .)
    myDmat <- lapply(1:length(dHmats), function(i){
      lapply(1:(q^2), function(j){
        eval_fn(dHmats[[i]][[j]], new_dataList[[i]], par_val, get_gradient = TRUE)
      }) %>% do.call(rbind, .)
    })

    traces <- sapply(1:q1, function(i){
      dH_val <- lapply(1:length(myDmat), function(j){
        matrix(myDmat[[j]][, i], q, q)
      }) %>% Reduce("+", .)
      matrix.trace(invH %*% (-dH_val))
    })
    fy1 <- -(dh_val - 0.5*traces)

    # derivative w.r.t. var-cov matrix parameters
    fy2 <- sapply(1:q0, function(i){
      dM_val <- evalMat(L2$dM[, , i], par_val = Lval)
      Dmat <- -mat%*%dM_val%*%mat # d invSIGMA/d xi
      XDmat <- diag(as.matrix(Bi_df[, -1])%*%Dmat%*%t(as.matrix(Bi_df[, -1])))
      if(distribution == "normal"){
        # from l_h
        gh1 <- sum(-0.5 * matrix.trace(mat%*% dM_val) - 0.5*XDmat)
        # from correction term log(|H|)
        DH <- Dmat*n
      } else if (distribution == "t-dist"){
        gh1 <- sum( 0.5* matrix.trace(invmat%*%Dmat) - (degree + q2)/(degree + deno)*0.5*XDmat)
        DH <- matrix(0, nrow=q2, ncol=q2)
        for(j in 1:n){
          XBi_df <- t(as.matrix(Bi_df[j, -1]))%*%as.matrix(Bi_df[j, -1])
          fg1 <- (degree + q2)*mat*(degree + deno[j]) - (degree + q2)*as.matrix(numo[j, ])%*%(2*numo[j, ]) # == numo1
          fg2 <- (degree+deno[j])^2
          Dg1 <- (degree*(degree + q2)+(degree + q2)*deno[j])*Dmat + (degree + q2)*mat*XDmat[j]-
            2*(degree + q2)*Dmat%*%XBi_df%*%mat - 2*(degree + q2)*mat%*%XBi_df%*%Dmat
          Dg2 <- as.numeric(2*(degree + deno[j])*XDmat[j])
          DH_j <- -(Dg1*fg2 - fg1*Dg2)/(fg2^2)
          DH <- DH + DH_j
        }
      }
      DH2 <- Matrix::bdiag(DH, diag(rep(0, q - q2))) %>% as.matrix()
      gh2 <- - 0.5*matrix.trace(invH%*%DH2)
      -(gh1+gh2)
    })

    c(fy1, fy2)
  }

  if(cov_method == "cholesky"){
    lower <- c(lower, rep(-Inf, q0))
    upper <- c(upper, rep(Inf, q0))
  } else if(cov_method == "spherical"){
    lower <- c(lower, rep(0, q0))
    upper <- c(upper, rep(pi, q0))
  }
  Lval0 <- initiate_covb_par(L2, invSIGMA, cov_method, correlation)
  str_val <- c(param, Lval0) %>% unlist()

  # check_gradient(str_val, ff, gr)

  res <- optim_iterate(str_val, ff, gr, lower = lower, upper = upper,
                       check = as.numeric(Silent), gr_tol = 1,
                       # extra info
                       param = param, other_param = other_param,
                       Hmats = Hmats, new_dataList = new_dataList,
                       RespLog = RespLog, Bi_df = Bi_df,
                       distribution = distribution, degree = degree,
                       dhlike = dhlike,
                       q = q, q1 = q1, q2 = q2, q0 = q0, n = n, L2 = L2,
                       dHmats = dHmats)
  disp_param <- res$gamma[1:q1]

  if(distribution == 'normal'){
    invmat <- evalMat(as.list(L2$M), par_val = res$gamma[-(1:q1)])
  } else if (distribution == 't-dist'){
    invmat <- evalMat(as.list(L2$M), par_val = res$gamma[-(1:q1)])*(degree - 2)/degree
  }
  mat <- solve(invmat)

  list(disp_param = disp_param, invSIGMA = mat, Lval = res$gamma[-(1:q1)])
}

initiate_covb_par <- function(L2, invSIGMA, cov_method, correlation){
  size_b <- nrow(invSIGMA)

  if(cov_method == "cholesky"){
    L_lower <- t(chol(solve(invSIGMA)))
    if(is.null(correlation)){
      Lval0 <- c()
      for(i in 1:size_b){
        Lval0 <- c(Lval0, L_lower[i, 1:i])
      }
    } else if(correlation == "independent"){
      Lval0 <- diag(L_lower)
    }
  } else if(cov_method == "spherical"){
    Lval0 <- rep(pi/2, length(L2$Mpar))
  }
  Vassign(L2$Mpar, Lval0)
}
