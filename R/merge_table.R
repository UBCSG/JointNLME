merge_table <- function(..., duplicate_col = 1){
  tables <- list(...)
  q <- ncol(tables[[1]])
  base_tb <- tables[[1]][, duplicate_col]
  
  for(i in setdiff(1:q, duplicate_col)){
    for(j in 1:length(tables)){
      base_tb <- cbind(base_tb, tables[[j]][, i])
    }
  }
  var_names <- colnames(tables[[1]])
  colnames(base_tb) <- c(var_names[duplicate_col],
                         rep(var_names[-duplicate_col], each = length(tables)))
  base_tb
}
