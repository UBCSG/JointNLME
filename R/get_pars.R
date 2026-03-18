
fixeff_par <- function(RespLog, random_effects, fixed_param, memList){
  q <- length(random_effects)
  p <- length(fixed_param)
  var_names <- names(fixed_param)
  # Derive -H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars = random_effects)

  # Derive -dH/d theta, where theta are fixed parameters
  dH_theta <- lapply(Hmats, function(hmat){
    all_cases <- strsplit(names(hmat), ",") %>% sapply(function(x){
      paste0(sort(x), collapse = ",")
    })
    lapply(1:p, function(i){
      covered_cases <- c()
      res <- rep("", q^2)
      for(j in 1:(q^2)){
        # cat("case=", all_cases[j],'\n')
        if(all_cases[j] %in% covered_cases){
          res[j] <- res[which(covered_cases == all_cases[j])[1]]
          # cat("exists already \n")
        } else {
          covered_cases <- c(covered_cases, all_cases[j])
          res[j] <- Vderiv(hmat[j], var_names[i])
          # cat("calculating \n")
        }
      }
      res
    })
  })

  # Derive d log_h/d theta,
  dfs <- lapply(RespLog, function(logh){
    lapply(1:p, function(i){ Vderiv(logh, var_names[i]) })
  })

  # Derive d (d log_h/dy) / d theta
  df_dy_dtheta <- lapply(1:length(RespLog), function(j){
      Vderiv(RespLog[[j]], memList[[j]]$response) %>%
        Vderiv(., var_names)
    })

  # Derive d^2 mu / db db^T
  H_mu <- getHmat(lapply(memList, function(mem_i){mem_i$reg_equation}), random_effects)

  list(Hmats = Hmats,
       dH_theta = dH_theta,
       dfs = dfs,
       df_dy_dtheta = df_dy_dtheta,
       H_mu = H_mu)
}

disp_par <- function(RespLog, deriv_vars, disp_param){
  q <- length(deriv_vars)
  # derivative of log-h w.r.t. param
  dhlike <- lapply(RespLog, function(logh){deriv(formula(paste("~", logh)), names(disp_param))})
  # get -H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars = deriv_vars)
  # derivative of -H w.r.t. param
  dHmats <- lapply(Hmats, function(hmat){
    lapply(1:q^2, function(i){
      deriv(formula(paste("~", hmat[i])), names(disp_param))
    })
  })
  list(dhlike = dhlike,
       Hmats = Hmats,
       dHmats = dHmats)
}
