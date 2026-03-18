# Function to create the 'trans' variable
create_trans <- function(df) {
  df$trans <- NA  # Initialize 'trans' column with NA
  df <- df %>% 
    arrange(ID, t)
  # Loop through each ID
  # for (id in 1:length(unique(df$ID))) {
  for (id in 33:33) {
    # Subset dataframe for each ID
    df_id <- df[df$ID == id, ]
    
    # Find the last row index with censor = 0 in decay phase
    last_decay_index <- max(which(df_id$phase_obs == "decay"))
    
    # Find the first row index with censor = 0 in rebound phase
    first_rebound_index <- min(which(df_id$phase_obs == "rebound" & df_id$censor == 0))
    
    # Set 'trans' values based on conditions
    df_id$trans <- NA
    df_id$trans[1: last_decay_index] <- 0
    df_id$trans[(first_rebound_index): nrow(df_id)] <- 1
    
    # Update 'trans' values in original dataframe
    df[df$ID == id, "trans"] <- df_id$trans
  }
  
  return(df)
}
 