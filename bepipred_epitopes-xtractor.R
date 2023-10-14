# Clear environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)

# Get a list of CSV files matching the pattern
bepipred_files <- list.files(pattern = "summary_.*\\.csv")

# Create an empty dataframe to store the combined data
combined_df <- data.frame()

# Loop through each CSV file and combine the data
for (file in bepipred_files) {
  # Read the current CSV file into a dataframe
  current_data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Combine the current dataframe with the combined_df
  combined_df <- rbind(combined_df, current_data)
}

# Filter to remain with only epitope probability score of 0.5 and above and columns on protein ID
bepipred_output <- combined_df %>%
  select(Entry, Position, AminoAcid, EpitopeProbability) %>%
  filter(EpitopeProbability >= 0.5)




# Initialize an empty dataframe to store the results
result_df <- data.frame(Entry = character(),
                        Start = numeric(),
                        End = numeric(),
                        Epitope = character(),
                        Length = numeric(),
                        stringsAsFactors = FALSE)

# Iterate through unique entries
unique_entries <- unique(bepipred_output$Entry)

for (entry in unique_entries) {
  entry_data <- bepipred_output[bepipred_output$Entry == entry, ]
  consecutive_amino_acids <- NULL
  current_start <- NA
  
  for (i in 1:nrow(entry_data)) {
    if (is.na(current_start)) {
      current_start <- entry_data$Position[i]
    }
    
    consecutive_amino_acids <- c(consecutive_amino_acids, entry_data$AminoAcid[i])
    
    if (i == nrow(entry_data) || entry_data$Position[i + 1] - entry_data$Position[i] > 1) {
      epitope <- paste(consecutive_amino_acids, collapse = "")
      if (nchar(epitope) > 0) {
        result_df <- rbind(result_df, data.frame(Entry = entry,
                                                 Start = current_start,
                                                 End = entry_data$Position[i],
                                                 Epitope = epitope,
                                                 Length = length(consecutive_amino_acids),
                                                 stringsAsFactors = FALSE))
      }
      current_start <- NA
      consecutive_amino_acids <- NULL
    }
  }
}

# Print the dataframe
print(result_df)


# remove epitopes with length less than 10 
library(tidyverse)

# Predicted b-cell linear epitopes 
bepipred_epitopes <- result_df %>%
  filter(Length>=10)

bepipred_epitopes$epitope_id <-  paste(bepipred_epitopes$Entry,
                                       bepipred_epitopes$Start, 
                                       bepipred_epitopes$End, sep = ":")

# rename column four to sequence

colnames(bepipred_epitopes)[4] <- "sequence" 

# Filter to remain with epitope id and sequence
bepipred_epitopes <- bepipred_epitopes %>%
  select(epitope_id, sequence)

# write out the results
#write_tsv(epitopes, file = "predicted_epitopes_bepipred.tsv")



##############################################################################
# epidope output combined 
#############################################################################
epidope_output <-  read_tsv("predicted_epitopes_epidope.csv")

# assign  column name for length
epidope_output$length <- (epidope_output$end - epidope_output$start) + 1 

# Assign epitope id 
epidope_output$epitope_id <- paste(epidope_output$`#Gene_ID`, 
                                   epidope_output$start, 
                                   epidope_output$end, 
                                   sep = "_")

# filter out epitopes ranging between 10-30
epidope_output <- epidope_output %>%
  select(epitope_id, sequence)

#write_tsv(epidope_output, file = "epidope_output.tsv") 


# combine epidope and bepipred_output 
combined_epitopes <- rbind(bepipred_epitopes, epidope_output)

# write_tsv 
write_tsv(combined_epitopes, file = 'combined_epitopes.tsv')


#################################################################################################################