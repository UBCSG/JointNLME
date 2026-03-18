#' Derive H matrix in profile h-likelihood functions
#'
#' Density of random effects is excluded here and will
#' be included separately while estimating parameters.
#' Note that this function returns -H.
#'
getHmat <- function(RespLog, pars){
  p <- length(pars)
  lapply(RespLog, function(logh){
    getHessian(logh, pars)
  })
}

#' Assign names to a vector
#'
Vassign <- function(name, value){
  dat <- data.frame(as.list(value))
  names(dat) <- name
  return(dat)
}

#' Calculate Hessian matrix
#'
#' @importFrom Deriv Simplify
getHessian <- function(lik1, pars){
  lik <- parse(text = lik1)
  q <- length(pars)
  result <- as.list(rep(NA, q^2))
  for(i in 1:q){
    for(j in 1:q){
      k <- (i - 1)*q + j
      result[[k]] <- Deriv::Simplify(D(D(lik, pars[i]), pars[j]))
      names(result)[k] <- paste(pars[i], pars[j], sep = ",")
    }
  }
  return(result)
}

#' This function returns a matrix with elements in string.
#' Moreover, it returns L(l)'L(l) = SIGMA, which is a spherical
#' parameterization with diagonal elements equal to 1.
#'
#' @importFrom Deriv Simplify
strMat <- function(q2, method = "cholesky"){
  L <- c()
  Mpar <- c()

  if(method == "spherical"){
    L <- rbind(L, c(1, rep(0, q2 - 1)))
    for(i in 2:q2){
      l0 <- 1
      Li <- c()
      for(j in 1:i){
        m0 <- paste('L', i, 2:min((j+1), i), sep = "")
        Mpar <- c(Mpar, m0)

        if(j < i){
          sin0 <- rep(c("sin(", "cos("), c(length(m0) - 1, 1))
        } else {
          sin0 <- rep("sin(",  length(m0))
        }
        l1 <- paste(sin0, m0, rep(")", length(m0)), sep = "")
        Li <- c(Li, paste(c(l0, l1), collapse = "*"))
      }
      L <- rbind(L, c(Li, rep(0, q2 - i)))
    }
  } else if(method == "cholesky"){
    for(i in 1:q2){
      m0 <- paste0("L", i, 1:i)
      L <- rbind(L, c(m0, rep(0, q2-i)))
      Mpar <- c(Mpar, m0)
    }
  }

    L <- t(L)

    M <- matrix(NA, q2, q2)  # M=L'L
    for(i in 1:q2){
      for(j in 1:q2){
        M[i, j] <- Deriv::Simplify(paste(L[, i], L[, j], sep = "*", collapse = "+"))
      }
    }

    if(method == "spherical"){ diag(M) <- "1" }

    Mpar <- unique(Mpar)
    Mp <- length(Mpar)
    dM <- array(NA, dim = c(q2, q2, Mp))
    for(i in 1:Mp){
      ttz <- unlist(lapply(M, function(x){ Vderiv(x, Mpar[i]) }))
      dM[,,i] <- matrix(as.character(ttz), q2, q2)
    }

    dL <- array(NA, dim = c(q2, q2, Mp))
    for(i in 1:Mp){
      ttz <- unlist(lapply(t(L), function(x){ Vderiv(x, Mpar[i]) }))
      dL[,,i] <- matrix(as.character(ttz), q2, q2)
    }

  return(list(M=M, Mpar=Mpar, dM=dM, L=t(L), dL=dL))
}

#'
#' @importFrom Deriv Simplify
strMat_ind <- function(q2, method = "cholesky"){
  L <- c()
  Mpar <- c()
  if(method == "cholesky"){
    for(i in 1:q2){
      m0 <- paste0("L", i, i)
      L <- rbind(L, c(rep(0, i-1), m0, rep(0, q2-i)))
      Mpar <- c(Mpar, m0)
    }
  }
  L <- t(L)

  M <- matrix(NA, q2, q2)  # M=L'L
  for(i in 1:q2){
    for(j in 1:q2){
      M[i, j] <- Deriv::Simplify(paste(L[, i], L[, j], sep = "*", collapse = "+"))
    }
  }

  Mpar <- unique(Mpar)
  Mp <- length(Mpar)
  dM <- array(NA, dim = c(q2, q2, Mp))
  for(i in 1:Mp){
    ttz <- unlist(lapply(M, function(x){ Vderiv(x, Mpar[i]) }))
    dM[,,i] <- matrix(as.character(ttz), q2, q2)
  }

  dL <- array(NA, dim = c(q2, q2, Mp))
  for(i in 1:Mp){
    ttz <- unlist(lapply(t(L), function(x){ Vderiv(x, Mpar[i]) }))
    dL[,,i] <- matrix(as.character(ttz), q2, q2)
  }

  return(list(M=M, Mpar=Mpar, dM=dM, L=t(L), dL=dL))
}

# Function to get decay and rebound phase only
get_decay_rebound <- function(df) {
  vl <- df$y_decay
  time <- df$time_decay
  nadir_idx <- which.min(vl)
  
  # Ensure there are enough points after nadir
  if (nadir_idx >= length(vl) - 3) return(df[1:nadir_idx, ])
  
  rebound_start_idx <- nadir_idx + which.max(diff(vl[(nadir_idx + 1):length(vl)]))
  
  diff_vl <- diff(vl)
  rolling_sd <- zoo::rollapply(diff_vl[rebound_start_idx:length(vl)], width = 3, sd, fill = NA)
  
  plateau_relative_idx <- which(rolling_sd < 0.1)[1]
  
  if (is.na(plateau_relative_idx)) {
    cutoff_idx <- length(vl)
  } else {
    cutoff_idx <- rebound_start_idx + plateau_relative_idx
  }
  
  df[1:cutoff_idx, ]
}


create_covb_param <- function(q, cov_method, correlation){
  if(is.null(correlation)){
    L2  <- strMat(q, method = cov_method)
  } else if(correlation == "independent"){
    if(cov_method == "cholesky"){
      L2  <- strMat_ind(q, method = cov_method)
    } else if(cov_method == "spherical"){
      L2 <- NULL # No need to estimate correlation matrix
    }
  }
  L2
}

#' Compute derivatives of an expression w.r.t
#' a vector of parameters respectively
#'
Vderiv <- function(lik1, pars){
  lik <- parse(text = lik1)
  q <- length(pars)
  result <- as.list(rep(NA, q))

  for(i in 1:q){
    result[[i]] <- Deriv::Simplify( D(lik, pars[i]) )
    names(result)[i] <- pars[i]
  }
  return(result)
}


#' Evaluate a matrix given texts
#'
#' @param mat a matrix of texts
#' @param data data used for evaluation
#' @param par_val parameter values needed for evaluation
#' @param byrow logical. If TRUE (the default) the matrix is filled by rows.
evalMat <- function(mat, data = NULL, par_val = NULL, byrow = TRUE){
  q <- sqrt(length(mat))
  if(q%%1 != 0){
    stop("mat must be a square matrix.")
  }
  par_val <- as.list(par_val)

  res <- lapply(mat, function(x){
    if(class(x) == "character"){
      res <- eval(parse(text = x), c(data, par_val)) %>% as.numeric()
    } else {
      res <- eval(x, c(data, par_val)) %>% as.numeric()
    }
    if(!is.null(data)){
      if(nrow(data) > 1 & length(res) == 1){
        # result is independent from data values
        res <- res*nrow(data)
      }
    }
    sum(res)
  })

  matrix(as.numeric(res), q, q, byrow = byrow)
}

eval_fn <- function(fn, data, par_val, get_gradient = TRUE){
  
  if(get_gradient){
    res <- attr(eval(fn, c(data, par_val)), "gradient")
    n <- nrow(res)
  } else {
    if(class(fn) %in% c("character", "list")){
      res <- eval(parse(text = fn), c(data, par_val))
    } else {
      res <- eval(fn, c(data, par_val))
    }
    n <- length(res)
  }

  if(!is.null(data)){
    if(n == 1 & nrow(data) > 1){
      # result is independent from data values
      res <- res*nrow(data)
    }
  }

  if(get_gradient){
    apply(res, 2, sum)
  } else {
    sum(res)
  }
}
# eval_fn <- function(fn, data, par_val, get_gradient = TRUE){
#   
#   if(get_gradient){
#     res <- attr(eval(fn, c(data, par_val)), "gradient")
#     n <- nrow(res)
#   } else {
#     if(class(fn) %in% c("character", "list")){
#       res <- eval(parse(text = fn), c(data, par_val))
#     } else {
#       res <- eval(fn, c(data, par_val))
#     }
#     n <- length(res)
#   }
# 
#   if(!is.null(data)){
#     if(n == 1 & nrow(data) > 1){
#       # result is independent from data values
#       res <- res*nrow(data)
#     }
#   }
# 
#   if(get_gradient){
#     apply(res, 2, sum)
#   } else {
#     sum(res)
#   }
# }

evalMat_row <- function(mat, data = NULL, par_val = NULL){
  par_val <- as.list(par_val)

  res <- lapply(mat, function(x){
    if(class(x) == "character"){
      res <- eval(parse(text = x), c(data, par_val)) %>% as.numeric()
    } else {
      res <- eval(x, c(data, par_val)) %>% as.numeric()
    }
    if(!is.null(data)){
      if(nrow(data) > 1 & length(res) == 1){
        # result is independent from data values
        res <- rep(res, nrow(data))
      }
    }
    res
  }) %>% do.call("cbind", .)

  res
}

eval_fn_row <- function(fn, data, par_val, get_gradient = TRUE){

  if(get_gradient){
    res <- attr(eval(fn, c(data, par_val)), "gradient")
    n <- nrow(res)
  } else {
    if(class(fn) %in% c("character", "list")){
      res <- eval(parse(text = fn), c(data, par_val))
    } else {
      res <- eval(fn, c(data, par_val))
    }
    n <- length(res)
  }

  if(!is.null(data)){
    if(n == 1 & nrow(data) > 1){
      # result is independent from data values
      res <- do.call("rbind", replicate(nrow(data), res, simplify = FALSE))
    }
  }
  res
}

check_gradient <- function(par, ff, gr, ..., delta  = 1e-4){
  f_val <- ff(par, ...)
  if(is.na(f_val)){
    stop("NaNs produced by ff().")
  }
  cat("ff = ", round(as.numeric(f_val), 3), "\n")
  g_val <- gr(par, ...)
  # cat("gr = ", round(as.numeric(g_val), 3), "\n")

  yy <- rep(NA, length(par))
  for(i in 1:length(par)){
    xx1 <- xx2 <- par
    xx1[i] <- xx1[i] + delta
    xx2[i] <- xx2[i] - delta
    yy[i] <- (ff(xx1, ...) - ff(xx2, ...))/(delta*2)
  }
  # cat("numerical gr =", round(as.numeric(yy), 3), "\n")
  # cat("gr diff = ", round(as.numeric(g_val - yy), 3), "\n")
  df <- data.frame(gr = g_val, num_gr = yy,
                   diff = g_val - yy)
  df <- round(df, 6)
  df$par <- names(par)
  df
}

num_gradient <- function(par, ff, ..., delta  = 1e-4){
  f_val <- ff(par, ...)
  if(is.na(f_val)){
    stop("NaNs produced by ff().")
  }
  cat("ff = ", round(as.numeric(f_val), 3), "\n")

  yy <- rep(NA, length(par))
  for(i in 1:length(par)){
    xx1 <- xx2 <- par
    xx1[i] <- xx1[i] + delta
    xx2[i] <- xx2[i] - delta
    yy[i] <- (ff(xx1, ...) - ff(xx2, ...))/(delta*2)
  }
  yy
}

prod_matrix3 <- function(A, B, C){
  A%*%B%*%C
}

prod_matrix2 <- function(A, B){
  A%*%B
}

neg_matrix <- function(A){
  -A
}

filter_data_by_id <- function(dat, subject_id, id){
  subset(dat, dat[, subject_id] == id)
}

rescale <- function(x, from = c(a, b), to = c(c, d)) {
  (x - from[1]) * (to[2] - to[1]) / (from[2] - from[1]) + to[1]
}
