predict <- function(object, ...){
  UseMethod("predict", object)
}

predict.hhjm <- function(jm_res){
  new_dataList <- lapply(jm_res$dataList, function(dat){dplyr::left_join(dat, jm_res$Bi, by = jm_res$subject_id)})
  y_hat <- list()
  for(i in 1:length(jm_res$memList)){
    y_hat[[jm_res$memList[[i]]$response]] <- with(new_dataList[[i]], with(as.list(c(jm_res$fixed_est, jm_res$disp_est)),
                                                         eval(parse(text = jm_res$memList[[i]]$reg_equation))))
  }
  y_hat
}
