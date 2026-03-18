simtable_withboot <- function(par_table, parsd_table, boot_sd, par_names, true_value, disppar_table, disppar_names, disp_value){
  nsim <- nrow(par_table)
  table <- data.frame(matrix(NA, nrow = length(true_value), ncol = 10))
  colnames(table) <- c("Parameter", "True Value", "Estimates", "SE_M", "SE_S", "SE_B", "Bias(%)", "rMSE(%)", "Coverage", "Coverage_B")
  table$Parameter <- par_names
  table$`True Value` <- true_value
  table$Estimates <- sapply(par_table, mean, na.rm = TRUE)
  table$SE_M <- sapply(parsd_table, mean, na.rm = TRUE)
  table$SE_S <- sapply(par_table, sd, na.rm = TRUE)
  table$SE_B <- sapply(boot_sd, mean, na.rm = TRUE)
  table$`Bias(%)` <- (table$Estimates - true_value)/true_value * 100
  bias_table <- sweep(par_table, 2, true_value)
  table$`rMSE(%)` <- apply(bias_table, 2, function(x){
    sqrt(mean(x^2))
  })/true_value * 100
  lower_CI <- par_table - 1.96 * parsd_table
  upper_CI <- par_table + 1.96 * parsd_table
  lower_CI <- sweep(lower_CI, 2, true_value, `<=`)
  upper_CI <- sweep(upper_CI, 2, true_value, `>=`)
  CI <- as.data.frame(lower_CI & upper_CI)
  table$Coverage <- sapply(CI, mean, na.rm = TRUE)
  lower_CI_B <- par_table - 1.96 * boot_sd
  upper_CI_B <- par_table + 1.96 * boot_sd
  lower_CI_B <- sweep(lower_CI_B, 2, true_value, `<=`)
  upper_CI_B <- sweep(upper_CI_B, 2, true_value, `>=`)
  CI_B <- as.data.frame(lower_CI_B & upper_CI_B)
  table$Coverage_B <- sapply(CI_B, mean, na.rm = TRUE)
  
  disptable <- data.frame(matrix(NA, nrow = length(disp_value), ncol = 10))
  colnames(disptable) <- c("Parameter", "True Value", "Estimates", "SE_M", "SE_S", "SE_B", "Bias(%)", "rMSE(%)", "Coverage", "Coverage_B")
  disptable$Parameter <- disppar_names
  disptable$`True Value` <- disp_value
  disptable$Estimates <- sapply(disppar_table, mean, na.rm = TRUE)
  disptable$SE_S <- sapply(disppar_table, sd, na.rm = TRUE)
  disptable$`Bias(%)` <- (disptable$Estimates - disp_value)/disp_value * 100
  dispbias_table <- sweep(disppar_table, 2, disp_value)
  disptable$`rMSE(%)` <- apply(dispbias_table, 2, function(x){
    sqrt(mean(x^2))
  })/disp_value * 100
  
  
  summary_table <- rbind(table, disptable)
  round_df(summary_table, 4)
}


round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
