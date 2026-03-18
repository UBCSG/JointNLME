get_sd_hlike <- function(fixed_param, random_effects,
                         other_param,
                         RespLog, new_dataList,
                         Bi_df, invSIGMA, T_end, last_decayobs,
                         disp_pars, 
                         Hmats3,
                         all_glmeObjects,
                         distribution = "normal", degree = 0){
  
  subject_id <- names(Bi_df)[1]
  random_effects <- names(Bi_df)[-1]
  q <- length(random_effects)
  n <- nrow(Bi_df)
  p <- length(fixed_param)
  uniqueID <- Bi_df[ ,1]
  # Hmats <- disp_pars$Hmats
  par_val <- c(fixed_param, other_param)
  
  # Derive -H matrix, where H is defined in the adjusted profile h-likelihood
  
  # evaluate the -H matrix for the random effects part
  if(distribution == "normal"){
    negH3 <- lapply(1:n, function(i){
      -invSIGMA %>% Matrix::bdiag(., diag(0, p+1, p+1)) %>% as.matrix()
    })
  } else if (distribution == "t-dist"){
    deno <- diag(as.matrix(Bi_df[, -1])%*%invSIGMA%*%t(Bi_df[, -1]))
    numo <- as.matrix(Bi_df[, -1])%*%invSIGMA
    negH3 <- lapply(1:n, function(i){
      numo1 <- (degree + q)*invSIGMA*(degree + deno[i]) -
        (degree + q)*as.matrix(numo[i, ])%*%(2*numo[i, ])
      -numo1/((degree + deno[i])^2) %>% Matrix::bdiag(., diag(0, p, p)) %>% as.matrix()
    })
  }
  
  # To approximate -H by E(-H|b) for NLMEs
  # Ref: Maengseok Noh and Youngjo Lee, 2008, Hierarchical-likelihood approach for nonlinear mixed-effects models
  exp_dataList <- new_dataList
  # for(i in 1:length(new_dataList)){
  #   for(j in 1:length(all_glmeObjects[[i]])){
  #     exp_dataList[[i]][all_glmeObjects[[i]][[j]]$response] <-
  #       eval_fn_row(all_glmeObjects[[i]][[j]]$reg_equation, new_dataList[[i]], par_val, get_gradient = F)
  #   }
  # }
  # hhh <- lapply(1:length(Hmats), function(i){
  #   evalMat_row(Hmats[[i]], exp_dataList[[i]], par_val) %>% as.data.frame() %>%
  #     cbind(exp_dataList[[i]][subject_id], .)
  # }) %>% do.call("rbind", .) %>%
  #   aggregate( as.formula(paste(". ~", subject_id)), ., sum)
  # 
  # 
  # B_mat <- NULL # matrix(0, n*q, n*q)
  # C_mat <- c() # matrix(0, n*q, p)
  # D_mat <- matrix(0, p, p)
  # for(k in 1:n){
  #   H_mat <- -(matrix(as.numeric(hhh[k, -1]), p+q, p+q) + negH3[[k]])
  #   if(is.null(B_mat)){
  #     B_mat <- H_mat[1:q, 1:q]
  #   } else {
  #     B_mat <- bdiag(B_mat, H_mat[1:q, 1:q])
  #   }
  #   C_mat <- rbind(C_mat, H_mat[1:q, (q+1):(q+p)])
  #   D_mat <- D_mat + H_mat[(q+1):(p+q), (q+1):(p+q)]
  # }
  
  # evaluate the H matrix
  hhh <- lapply(1:length(Hmats3), function(i){
    evalMat_row(Hmats3[[i]], exp_dataList[[i]], par_val) %>% as.data.frame() %>%
      cbind(exp_dataList[[i]][subject_id], .)
  }) %>% do.call("rbind", .) %>%
    aggregate( as.formula(paste(". ~", subject_id)), ., sum)
  
  B_mat <- NULL # matrix(0, n*q, n*q)
  C_mat <- c() # matrix(0, n*q, p)
  D_mat <- matrix(0, p+1, p+1)
  for(k in 1:n){
    H_mat <- -(matrix(as.numeric(hhh[k, -1]), p+q+1, p+q+1) + negH3[[k]])
    if(is.null(B_mat)){
      B_mat <- H_mat[1:q, 1:q]
    } else {
      B_mat <- bdiag(B_mat, H_mat[1:q, 1:q])
    }
    C_mat <- rbind(C_mat, H_mat[1:q, (q+1):(q+p+1)])
    D_mat <- D_mat + H_mat[(q+1):(p+q+1), (q+1):(p+q+1)]
  }
  
  Hval <- rbind(cbind(B_mat, C_mat), cbind(t(C_mat), D_mat))
  sd2 <- solve(Hval) %>% diag() %>% tail(p)
  sqrt(sd2)
}
