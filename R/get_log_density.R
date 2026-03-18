#' Get log density function of (nonlinear) (generalized) linear model or survial model
#'
get_log_density <- function(model_object){
  if(model_object$modeltype == "nlme"){
    get_log_density_glm(model_object)
  } else if (model_object$modeltype == "survival"){
    get_log_density_ph(model_object)
  }
}

#' #' Get joint log density function
#' get_joint_density <- function(memList){
#'   if(!is.null(memList)){
#'     long_config <- lapply(memList, get_log_density_glm)
#'     long_log_hlike <- lapply(long_config, function(x){x$log_density}) %>% paste(collapse = "+")
#'   } else {
#'     long_config <- long_log_hlike <- NULL
#'   }
#'   long_log_hlike
#' }

#' Get log density function of (generalized) linear model
#'
get_log_density_glm <- function(model_object){
  if(is.null(model_object)){
    return(NULL)
  }

  additional_param <- NULL
  left_censored <- !is.null(model_object$left_censoring)

  if(!left_censored){

    if(model_object$distribution == "binomial"){

      res <- paste(model_object$response, "*(",
                   model_object$reg_equation, ")-", "log(1+exp(",
                   model_object$reg_equation, "))")

    } else if(model_object$distribution == "normal"){

      sigma <- paste0(model_object$response, "_sigma")
      res <- paste("- 0.5*(", model_object$response, "- (",
                   model_object$reg_equation, "))^2/",
                   sigma, "^2-log(", sigma, ")-0.5*log(2*pi)")
      additional_param <- c(additional_param, sigma)

    } else if(model_object$distribution == "poisson"){

      res <- paste(model_object$response, "*(",
                   model_object$reg_equation, ")-exp(",
                   model_object$reg_equation, ")-log(factorial(",
                   model_object$response, "))")
    }

  } else {

    sigma <- paste0(model_object$response, "_sigma")
    log_density <- paste("- 0.5*(", model_object$response,
                         "- (", model_object$reg_equation, "))^2/",
                         sigma, "^2-log(", sigma, ")-0.5*log(2*pi)")
    Cresp <- model_object$left_censoring$indicator

    if(model_object$left_censoring$method == "tobit"){
      # limit_val <- model_object$response
      limit_val <- model_object$left_censoring$limit_value
    } else if(model_object$left_censoring$method == "truncated"){
      stopifnot(!is.null(model_object$left_censoring$limit_value))
      limit_val <- model_object$left_censoring$limit_value
    }

    CstdNmu <- paste("(", limit_val, "-(",
                     model_object$reg_equation, "))/", sigma)
    leftint <- paste("1/(1+exp(-1.702*", CstdNmu, "))")

    if(model_object$left_censoring$method == "tobit"){
      # A Tobit model, assuming left-censored data follow normal distribution
      res <- paste( "(", log_density, ")*(1-", Cresp, ") +log(", leftint, ")*",Cresp )

    } else if(model_object$left_censoring$method == "truncated"){
      # Model uncensored data only, assumed following truncted normal distribution
      res <- paste("(", log_density, "-log(1-", leftint, "))*(1-", Cresp, ")")
    }

    additional_param <- c(additional_param, sigma)
  }

  list(log_density = paste("(", res, ")"),
       additional_param = additional_param)
}

#' Get log density function of propotional hazard model
#'
get_log_density_ph <- function(model_object){
  additional_param <- NULL
  if(is.null(model_object$distribution)){
    log_density <- paste(model_object$event, "*( log(hazard0) +",
                         model_object$reg_equation, ") - exp(",
                         model_object$reg_equation, ")*cum_hazard0")

    additional_param <- c("hazard0", "cum_hazard0")
    
  } else if(model_object$distribution == "weibull"){
    # Weibull PH model
    # log_density <- paste(model_object$event, "*( Wlogscale+log(Wshape)+log(",
    #                      model_object$response, ")*(Wshape-1) +",
    #                      model_object$reg_equation, ") - exp(Wlogscale+",
    #                      model_object$reg_equation, ")*",
    #                      model_object$response, "^Wshape")
    # log_density <- paste(model_object$event, "*( log(Wscale)+log(Wshape)+log(",
    #                      model_object$response, ")*(Wshape-1) +",
    #                      model_object$reg_equation, ") - exp(log(Wscale) + ",
    #                      model_object$reg_equation, ")*",
    #                      model_object$response, "^Wshape")
    log_density <- paste(model_object$event, "*( log(Wscale)+log(Wshape)+log(",
                         model_object$response, ")*(Wshape-1) +",
                         model_object$reg_equation, ") - exp(",
                         model_object$reg_equation, ") * Wscale *",
                         model_object$response, "^Wshape")
    # Wscale = sigma = root(1/lambda, gamma)
    # Wshape = a = gamma
    # additional_param <- c("Wscale", "Wshape")
  } else if (model_object$distribution == "log-normal"){
    sigma <- paste0(model_object$response, "_sigma")
    upper_bound <- paste("( log(", model_object$truncated$upper, "- ", model_object$truncated$lower, ") - (",
                         model_object$reg_equation, "))/", sigma)
    upperint <- paste("1/(1+exp(-1.702*", upper_bound, "))")
    log_density <- paste("-log(", model_object$response,") - 0.5*(log(", model_object$response, ") - (",
                         model_object$reg_equation, "))^2/",
                         sigma, "^2-log(", sigma, ")-0.5*log(2*pi) - log(", upperint, ")")
    #   upper_bound <- paste("( log(", model_object$truncated$upper, ") - (",
    #                        model_object$reg_equation, "))/", sigma)
    #   upperint <- paste("1/(1+exp(-1.702*", upper_bound, "))")
    #   lower_bound <- paste("( log(", model_object$truncated$lower, ") - (",
    #                        model_object$reg_equation, "))/", sigma)
    #   lowerint <- paste("1/(1+exp(-1.702*", lower_bound, "))")
    # if (model_object$truncated$lower == 0){
    #   upper_bound <- paste("( log(", model_object$truncated$upper, "- ", model_object$truncated$lower, ") - (",
    #                        model_object$reg_equation, "))/", sigma)
    #   upperint <- paste("1/(1+exp(-1.702*", upper_bound, "))")
    #   log_density <- paste("-log(", model_object$response,") - 0.5*(log(", model_object$response, ") - (",
    #                        model_object$reg_equation, "))^2/",
    #                        sigma, "^2-log(", sigma, ")-0.5*log(2*pi) - log(", upperint, ")")
    # } else{
    #   log_density <- paste("-log(", model_object$response,") - 0.5*(log(", model_object$response, ") - (",
    #                        model_object$reg_equation, "))^2/",
    #                        sigma, "^2-log(", sigma, ")-0.5*log(2*pi) - log(", upperint ," - ", lowerint, ")")
    # }

    additional_param <- c(additional_param, sigma)
  }

  list(log_density = log_density,
       additional_param = additional_param)
}

