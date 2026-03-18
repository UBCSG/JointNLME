get_sd <- function(fixed_param, random_effects,
                   other_param,
                   RespLog, new_dataList,
                   Bi_df, invSIGMA, T_end, last_decayobs,
                   Hmats2,
                   all_glmeObjects,
                   distribution = "normal", degree = 0){
  
  subject_id <- names(Bi_df)[1]
  random_effects <- names(Bi_df)[-1]
  q <- length(random_effects)
  n <- nrow(Bi_df)
  p <- length(fixed_param)
  uniqueID <- Bi_df[ ,1]
  par_val <- c(fixed_param, other_param)

  # evaluate the H matrix
  hhh <- lapply(1:length(Hmats2), function(i){
    evalMat_row(Hmats2[[i]], new_dataList[[i]], par_val) %>% as.data.frame() %>%
      cbind(new_dataList[[i]][subject_id], .)
  }) %>% do.call("rbind", .) %>%
    aggregate( as.formula(paste(". ~", subject_id)), ., sum)

  H_mat_final <- matrix(0, p, p)
  for (z in 1:n){
    H_mat <- -(matrix(as.numeric(hhh[z, -1]), p, p))
    H_mat_final <- H_mat_final + H_mat
  }
  
  sd2 <- solve(H_mat_final) %>% diag()
  sqrt(sd2)
}
