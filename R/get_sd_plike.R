# eval_H_b <- function(RespLog, random_effects, dataList,
#                      Bi_df, invSIGMA, par_val,
#                      distribution, degree){
#   n <- nrow(Bi_df)
#   q <- length(random_effects)
#
#   # Derive H matrix, where H is defined in the adjusted profile h-likelihood
#   Hmats <- getHmat(RespLog, pars = random_effects)
#   deno <- diag(as.matrix(Bi_df[, -1])%*%invSIGMA%*%t(Bi_df[, -1]))
#   numo <- as.matrix(Bi_df[, -1])%*%invSIGMA
#   # evaluate the -H matrix for the random effects part
#   if(distribution == "normal"){
#     negH3 <- -invSIGMA*n
#   } else if (distribution == "t-dist"){
#     negH3 <- lapply(1:n, function(i){
#       numo1 <- (degree + q)*invSIGMA*(degree + deno[i]) -
#         (degree + q)*as.matrix(numo[i, ])%*%(2*numo[i, ])
#       -numo1/((degree + deno[i])^2)
#     }) %>% Reduce('+', .)
#   }
#   negH <- lapply(1:length(Hmats), function(i){
#     evalMat(Hmats[[i]], dataList[[i]], par_val = par_val)
#   }) %>% Reduce("+", .)
#   -(negH + negH3)
# }
#
# eval_H_theta <- function(RespLog, fixed_param, other_param, dataList){
#   par_val <- c(fixed_param, other_param)
#   Hmats <- getHmat(RespLog, pars = names(fixed_param))
#   negH <- lapply(1:length(Hmats), function(i){
#     evalMat(Hmats[[i]], dataList[[i]], par_val = par_val)
#   }) %>% Reduce("+", .)
#   -negH
# }
#
# eval_dHb_dtheta <- function(RespLog, random_effects,
#                             fixed_param, dataList, par_val){
#   q <- length(random_effects)
#   p <- length(fixed_param)
#
#   Hmats <- getHmat(RespLog, pars = random_effects)
#   dH_i <- array(0, dim = c(q, q, p))
#   for(i in 1:p){
#     dH_i[, , i] <- lapply(1:length(Hmats), function(k){
#       lapply(Hmats[[k]], function(x){
#         D(x, names(fixed_param)[i])
#       }) %>%
#         evalMat(dataList[[k]], par_val = par_val)
#     }) %>% Reduce("+", .)
#   }
#   -dH_i
# }
#
# eval_d2Hb_d2theta <- function(RespLog, random_effects,
#                             fixed_param, dataList, par_val){
#   q <- length(random_effects)
#   p <- length(fixed_param)
#   Hmats <- getHmat(RespLog, pars = random_effects)
#   dH_ij <- array(0, dim = c(q, q, p, p))
#   for(i in 1:p){
#     for(j in i:p){
#       Dij <- lapply(1:length(Hmats), function(k){
#         Vderiv(Hmats[[k]], names(fixed_param)[i]) %>%
#           Vderiv(names(fixed_param)[j]) %>%
#           evalMat(dataList[[k]], par_val = par_val)
#       }) %>% Reduce("+", .)
#       dH_ij[,,i,j] <- dH_ij[,,j,i] <- Dij
#     }
#   }
#   -dH_ij
# }
#
# get_sd_plike <- function(fixed_param, random_effects,
#                    other_param,
#                    RespLog, dataList,
#                    Bi_df, invSIGMA,
#                    # lambda = 0,
#                    distribution = "normal", degree = 0){
#
#   subject_id <- names(Bi_df)[1]
#   random_effects <- names(Bi_df)[-1]
#   q <- length(random_effects)
#   n <- nrow(Bi_df)
#   p <- length(fixed_param)
#   dataList <- lapply(dataList, function(dat){dplyr::left_join(dat, Bi_df, by = subject_id)})
#   par_val <- c(fixed_param, other_param)
#
#   Hval_theta <- eval_H_theta(RespLog, fixed_param, other_param, dataList)
#
#   Hval_b <- eval_H_b(RespLog, random_effects, dataList,
#                      Bi_df, invSIGMA, par_val, distribution, degree)
#   invH <- solve(Hval_b)
#   dH_i <- eval_dHb_dtheta(RespLog, random_effects, fixed_param, dataList, par_val)
#   d2H_ij <- eval_d2Hb_d2theta(RespLog, random_effects, fixed_param, dataList, par_val)
#
#   Hessian1 <- -Hval_theta
#   Hessian2 <- matrix(NA, p, p)
#   for(i in 1:p){
#     for(j in i:p){
#       Hessian2[i, j] <- Hessian2[j, i] <-
#         matrix.trace(invH%*%d2H_ij[,,i,j] - invH%*%dH_i[,,j]%*%invH%*%dH_i[,,i])
#     }
#   }
#   Hessian <- Hessian1 - 0.5*Hessian2
#   covMat <- -solve(Hessian)
#   sqrt(diag(covMat))
# }
