filter_dataList <- function(dataList, subject_id, id){
  lapply(dataList, function(dat){subset(dat, dat[, subject_id] == id)})
}

cbind_dataList <- function(dataList, Bi){
  lapply(dataList, function(dat){cbind(dat, Bi)})
}

get_uniqueID <- function(dataList, subject_id){
  lapply(dataList, function(dat){dat[, subject_id]}) %>% unlist %>% unique()
}

get_gr_logh <- function(RespLog, random_effects){
  lapply(RespLog, function(logh){deriv(formula(paste("~", logh)), random_effects)}) 
}

name_df <- function(data, colname){
  data <- as.data.frame(data)
  names(data) <- colname
  data
}
