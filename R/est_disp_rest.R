#' Estimate the dispersion parameters in the joint models,
#' by maximizing the adjusted profile h-likelihood function.
#'
#' @importFrom Matrix bdiag
#' @importFrom dplyr select
#' @importFrom matrixcalc matrix.trace
#' @importFrom dplyr left_join
est_disp_rest <- function(param, other_param,
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

  uniqueID <- unique(Bi_df[, 1])

  if(distribution == 'normal'){
    nD3 <- lapply(1:n, function(k){
      Matrix::bdiag( -invSIGMA, diag(0, q-q2, q-q2))
    })
  } else if(distribution == 't-dist'){
    numo <- as.matrix(Bi_df[, -1]) %*% invSIGMA
    deno <- numo%*%t(Bi_df[, -1]) %>% diag()
    nD3 <- list()
    for(i in 1:n){
      numo1 <- (degree + q2)*invSIGMA*(degree + deno[i]) - (degree + q2)*as.matrix(numo[i, ])%*%(2*numo[i, ])
      nD3_i <- -numo1/((degree + deno[i])^2)
      nD3[[i]] <- Matrix::bdiag( nD3_i, diag(0, q-q2, q-q2))
    }
  }

  # Return the negative value of the adjusted profile h-likelihood
  ff <- function(xx, ...){
    par_val <- c(Vassign(names(param), xx), other_param)

    #component by random effect
    # evaluate the H matrix
    hhh <- lapply(1:length(Hmats), function(i){
      evalMat_row(Hmats[[i]], new_dataList[[i]], par_val) %>%
        cbind(new_dataList[[i]][subject_id], .)
    }) %>% do.call("rbind", .) %>%
      aggregate( as.formula(paste(". ~", subject_id)), ., sum)
    A <- diag(0)
    C <- c()
    D <- matrix(0, q-q2, q-q2)
    for(k in 1:n){
      nD <- matrix(as.numeric(hhh[k, -1]), q, q)
      H_i <- - (nD + as.matrix(nD3[[k]]))
      A <- bdiag(A, H_i[1:q2, 1:q2])
      C <- rbind(C, H_i[c(1:q2), -c(1:q2)])
      D <- D + H_i[-c(1:q2), -c(1:q2)]
    }
    A <- as.matrix(A)
    A2 <- D - t(C)%*%solve(A)%*%C
    adj_val <- -0.5*( log(det(A/(2*pi))) + log(det(A2/(2*pi))) )
    # H <- rbind(cbind(A, C) , cbind(t(C), D))
    # adj_val <- -0.5*log(det(H/(2*pi)))

    # evaluate the h-likelihood
    lhlike <- eval_log_hlike(RespLog, new_dataList, par_val, as.matrix(Bi_df[, -1]),
                             invSIGMA, distribution, degree)
    # evaluate the negative adjusted profile h-likelihood
    -(lhlike + adj_val)
  }

  gr <- function(xx, ...){
    par_val <- c(Vassign(names(param), xx), other_param)
    # evaluate the H matrix
    hhh <- lapply(1:length(Hmats), function(i){
      evalMat_row(Hmats[[i]], new_dataList[[i]], par_val) %>%
        cbind(new_dataList[[i]][subject_id], .)
    }) %>% do.call("rbind", .) %>%
      aggregate( as.formula(paste(". ~", subject_id)), ., sum)
    A <- diag(0)
    C <- c()
    D <- matrix(0, q-q2, q-q2)
    for(k in 1:n){
      nD <- matrix(as.numeric(hhh[k, -1]), q, q)
      H_i <-  - (nD + as.matrix(nD3[[k]]))
      A <- bdiag(A, H_i[1:q2, 1:q2])
      C <- rbind(C, H_i[c(1:q2), -c(1:q2)])
      D <- D + H_i[-c(1:q2), -c(1:q2)]
    }
    A <- as.matrix(A)
    H <- rbind(cbind(A, C) , cbind(t(C), D))
    invH <- solve(H) %>% as.matrix()
    # comp1 <- solve(D)%*%t(C)
    # mat_e1 <- solve(A - C%*%comp1)
    # mat_e2 <- -mat_e1%*%t(comp1)
    # mat_e4 <- solve(D) + comp1%*%mat_e1%*%t(comp1)
    # invH <- rbind(cbind(mat_e1, mat_e2), cbind(t(mat_e2), mat_e4)) %>% as.matrix()

    # derivative w.r.t. dispersion parameters
    # part 1:
    dh_val <- lapply(1:length(dhlike), function(i){
      eval_fn(dhlike[[i]], new_dataList[[i]], par_val, get_gradient = TRUE)
    }) %>% Reduce("+", .)
    # part 2:
    myDmat <- lapply(1:(q^2), function(j){ # for each entry of H
      # calculate d Hij / d par by subject
      ggg <- lapply(1:length(dHmats), function(i){
        hhh <- eval_fn_row(dHmats[[i]][[j]], new_dataList[[i]], par_val, get_gradient = TRUE)
        hhh <- cbind(new_dataList[[i]][subject_id], hhh)
      }) %>% do.call("rbind", .) %>%
        aggregate( as.formula(paste(". ~", subject_id)), ., sum)
      cbind(h_idx = j, ggg)
    }) %>% do.call("rbind", .)

    traces <- sapply(1:q1, function(i){
      A <- diag(0)
      C <- c()
      D <- matrix(0, q-q2, q-q2)
      for(k in 1:n){
        dH_val <- subset(myDmat, myDmat[, subject_id] == uniqueID[k]) %>%
          arrange(h_idx) %>% .[, names(param)[i]] %>%
          matrix(., q, q)
        A <- bdiag(A, dH_val[1:q2, 1:q2])
        C <- rbind(C, dH_val[c(1:q2), -c(1:q2)])
        D <- D + dH_val[-c(1:q2), -c(1:q2)]
      }
      dH_val_all <- rbind(cbind(A, C), cbind(t(C), D))
      matrix.trace(as.matrix(invH %*% (-dH_val_all)))
    })
    -(dh_val - 0.5*traces)
  }

  str_val <- unlist(param)
  # check_gradient(str_val, ff, gr) %>% print
  # xx <- rnorm(length(str_val), str_val, 0.1) %>% pmin(., upper) %>% pmax(., lower)
  # ttt <- check_gradient(xx, ff, gr, delta = 1e-10)

  res <- optim_iterate_nlminb(str_val, ff, gr, lower = lower, upper = upper,
                              check = as.numeric(Silent), gr_tol = 1,
                              # extra info
                              param = param, other_param = other_param,
                              Hmats = Hmats, new_dataList = new_dataList,
                              RespLog = RespLog, Bi_df = Bi_df,
                              distribution = distribution, degree = degree,
                              dhlike = dhlike,
                              q = q, q1 = q1, q2 = q2, q0 = q0, n = n, L2 = L2,
                              dHmats = dHmats,
                              uniqueID = uniqueID)
  res$gamma
}
