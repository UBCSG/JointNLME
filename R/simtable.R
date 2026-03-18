simtable <- function(par_table2, parsd_table2, par_names, true_value, disppar_table2, disppar_names, disp_value, boot_sd2 = NULL){
  nsim <- nrow(par_table2)
  if (is.null(boot_sd2)){
    table <- data.frame(matrix(NA, nrow = length(true_value), ncol = 8))
    disptable <- data.frame(matrix(NA, nrow = length(disp_value), ncol = 8))
    colnames(table) <- colnames(disptable) <-c("Parameter", "True Value", "Estimates", "SE_M", "SE_S", "Bias(%)", "rMSE(%)", "CR_M")
  } else{
    table <- data.frame(matrix(NA, nrow = length(true_value), ncol = 10))
    disptable <- data.frame(matrix(NA, nrow = length(disp_value), ncol = 10))
    colnames(table) <- colnames(disptable) <- c("Parameter", "True Value", "Estimates", "SE_M", "SE_B", "SE_S", "Bias(%)", "rMSE(%)", "CR_M", "CR_B")
    table$SE_B <- sapply(boot_sd2, mean, na.rm = TRUE)
  }
  table$Parameter <- par_names
  table$`True Value` <- true_value
  table$Estimates <- sapply(par_table2, mean, na.rm = TRUE)
  table$SE_M <- sapply(parsd_table2, mean, na.rm = TRUE)
  table$SE_S <- sapply(par_table2, sd, na.rm = TRUE)
  table$`Bias(%)` <- (table$Estimates - true_value)/true_value * 100
  bias_table <- sweep(par_table2, 2, true_value)
  table$`rMSE(%)` <- apply(bias_table, 2, function(x){
    sqrt(mean(x^2))
  })/abs(true_value) * 100
  lower_CI <- par_table2 - 1.96 * parsd_table2
  upper_CI <- par_table2 + 1.96 * parsd_table2
  lower_CI <- sweep(lower_CI, 2, true_value, `<=`)
  upper_CI <- sweep(upper_CI, 2, true_value, `>=`)
  CI <- as.data.frame(lower_CI & upper_CI)
  table$CR_M <- sapply(CI, mean, na.rm = TRUE)
  
  if (!is.null(boot_sd2)){
    lower_CI_B <- par_table2 - 1.96 * boot_sd2
    upper_CI_B <- par_table2 + 1.96 * boot_sd2
    lower_CI_B <- sweep(lower_CI_B, 2, true_value, `<=`)
    upper_CI_B <- sweep(upper_CI_B, 2, true_value, `>=`)
    CI_B <- as.data.frame(lower_CI_B & upper_CI_B)
    table$CR_B <- sapply(CI_B, mean, na.rm = TRUE)
  }

  disptable$Parameter <- disppar_names
  disptable$`True Value` <- disp_value
  disptable$Estimates <- sapply(disppar_table2, mean, na.rm = TRUE)
  disptable$SE_S <- sapply(disppar_table2, sd, na.rm = TRUE)
  disptable$`Bias(%)` <- (disptable$Estimates - disp_value)/disp_value * 100
  dispbias_table <- sweep(disppar_table2, 2, disp_value)
  disptable$`rMSE(%)` <- apply(dispbias_table, 2, function(x){
    sqrt(mean(x^2))
  })/abs(disp_value) * 100

  
  summary_table <- rbind(table, disptable)
  round_df(summary_table, 4)
}


round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
