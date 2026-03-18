#' Check if the inputs of memList are valid
#'
check_memList <- function(memList){
  responses <- lapply(memList,
                      function(x){x$response}) %>% unlist()
  duplicate_name <- any(table(responses) > 1)
  if(duplicate_name){
    stop("There are duplicated response names. Rename the variables and try again.")
  } else {
    print("Ready to go!")
  }
}
