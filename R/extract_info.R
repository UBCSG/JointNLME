#'
extract_random_effects <- function(all_glmeObjects){
  sapply(all_glmeObjects, function(x){x$random_effects}) %>% unlist %>%
    as.character() %>% unique
}

#'
extract_fixed_param <- function(all_glmeObjects){
  fixed_param <- lapply(all_glmeObjects, function(x){
    res <- Vassign(x$fixed_param$names, x$fixed_param$start_values)
  }) %>% unlist() # %>% do.call(c, .)
  fixed_param <- fixed_param[unique(names(fixed_param))]
  # fixed_lower <- rep(-Inf, length(fixed_param))
  # fixed_upper <- rep(Inf, length(fixed_param))
  fixed_lower <- lapply(all_glmeObjects, function(x){x$fixed_param$lower}) %>% unlist
  fixed_upper <- lapply(all_glmeObjects, function(x){x$fixed_param$upper}) %>% unlist
  names(fixed_lower) <- names(fixed_upper) <- names(fixed_param)

  list(fixed_param = fixed_param,
       lower = fixed_lower,
       upper = fixed_upper)
}

#'
extract_disp_param <- function(all_glmeObjects){
  # dispersion parameters from coefficients of random effects
  disp_param <- lapply(all_glmeObjects, function(x){
    Vassign(x$disp_param$names, x$disp_param$start_values)
  })
  disp_lower <- lapply(all_glmeObjects, function(x){x$disp_param$lower}) %>% unlist
  disp_upper <- lapply(all_glmeObjects, function(x){x$disp_param$upper}) %>% unlist

  # dispersion parameters from residual terms in NLME models
  # sigma_param <- lapply(all_glmeObjects, function(x){
  #   if(x$modeltype=="nlme"){get_log_density_glm(x)$additional_param}
  # }) %>% unlist() %>% Vassign(., lapply(all_glmeObjects, function(x){ x$sigma }) %>% unlist())

  sigma_param <- lapply(all_glmeObjects, function(x){
    get_log_density(x)$additional_param
  }) %>% unlist() %>% Vassign(., lapply(all_glmeObjects, function(x){ x$sigma }) %>% unlist())
  
  sigma_lower <- rep(0, length(sigma_param))
  sigma_upper <- rep(Inf, length(sigma_param))

  list(# disp_param = do.call(c, c(disp_param, sigma_param)),
    disp_param = c(disp_param, sigma_param) %>% unlist(),
       lower = c(disp_lower, sigma_lower),
       upper = c(disp_upper, sigma_upper))
}

#'
extract_trans_param <- function(trans){
  # dispersion parameters from coefficients of random effects
  fixed_param <- Vassign(trans$mem$fixed_param$names, trans$mem$fixed_param$start_values) %>% as.list()
  fixed_lower <- trans$mem$fixed_param$lower
  fixed_upper <- trans$mem$fixed_param$upper
  disp_param <- Vassign(trans$mem$disp_param$names, trans$mem$disp_param$start_values) %>% as.list()
  disp_lower <- trans$mem$disp_param$lower
  disp_upper <- trans$mem$disp_param$upper
  

  list(param = c(fixed_param, disp_param),
       lower = c(fixed_lower, disp_lower),
       upper = c(fixed_upper, disp_upper))
}

#'
initiate_invSIGMA <- function(all_glmeObjects, randeff_info){
  diag_elem <- lapply(all_glmeObjects, function(x){
    if(!is.null(x$random_param)){
      x$random_param$sd^2
    } else {
      rep(1, length(x$random_effects))
    }
  }) %>% unlist()
  q <- length(diag_elem)
  SIGMA <- diag(diag_elem, q, q)
  if(randeff_info$distribution == 't-dist'){
    SIGMA <- diag(1, q, q)*(randeff_info$degree - 2)/randeff_info$degree
  }
  solve(SIGMA)
}
