
# clear environment 
 rm(list=ls())

# load packages 
 
 library(readr)
 library(tidyverse)

 
 # load sample data  
 
 bepipred_output <- read_csv('summary_64ECCD8C0000485EBF6CF1D2.csv', show_col_types = FALSE)


 # Filter to remain with only epitope probability score of 0.5 and above and columns on protein ID
 
 bepipred_output <- bepipred_output %>%
   select(Entry, Position, AminoAcid, EpitopeProbability)%>%
   filter(EpitopeProbability>=0.5)

 ################################################################################################################
#Bepipred_loops
 ###############################################################################################################
 
 # Create an empty dataframe to store the results
 result_df <- data.frame(Entry = character(),
                         Start = numeric(),
                         End = numeric(),
                         Epitope = character(),
                         Length = numeric(),
                         stringsAsFactors = FALSE)
 
 current_consecutive_amino_acids <- NULL
 
 for (i in 1:nrow(bepipred_output)) {
   current_consecutive_amino_acids <- c(current_consecutive_amino_acids, bepipred_output[i, 3])
   
   # Check if the next amino acid should be included in the epitope
   if (i == nrow(bepipred_output) || bepipred_output[i+1, 2] - bepipred_output[i, 2] > 1) {
     epitope <- paste(current_consecutive_amino_acids, collapse = "")
     
     # Add the epitope information to the dataframe if it's not empty
     if (nchar(epitope) > 0) {
       result_df <- rbind(result_df,
                          data.frame(Entry = bepipred_output[i, 1],
                                     Start = bepipred_output[i - length(current_consecutive_amino_acids) + 1, 2],
                                     End = bepipred_output[i, 2],
                                     Epitope = epitope,
                                     Length = length(current_consecutive_amino_acids)))
     }
     
     current_consecutive_amino_acids <- NULL
   }
 }
 
 # Print the dataframe
 print(result_df)
 
 #######################################################################################################################
 # Bepipred function
 
######################################################################################################################
# make a function 
 
 bepipred.extract=function(bepipred_output){
   # Create an empty dataframe to store the results
   result_df <- data.frame(Entry = character(),
                           Start = numeric(),
                           End = numeric(),
                           Epitope = character(),
                           Length = numeric(),
                           stringsAsFactors = FALSE)
   
   current_consecutive_amino_acids <- NULL
   
   for (i in 1:nrow(bepipred_output)) {
     current_consecutive_amino_acids <- c(current_consecutive_amino_acids, bepipred_output[i, 3])
     
     # Check if the next amino acid should be included in the epitope
     if (i == nrow(bepipred_output) || bepipred_output[i+1, 2] - bepipred_output[i, 2] > 1) {
       epitope <- paste(current_consecutive_amino_acids, collapse = "")
       
       # Add the epitope information to the dataframe if it's not empty
       if (nchar(epitope) > 0) {
         result_df <- rbind(result_df,
                            data.frame(Entry = bepipred_output[i, 1],
                                       Start = bepipred_output[i - length(current_consecutive_amino_acids) + 1, 2],
                                       End = bepipred_output[i, 2],
                                       Epitope = epitope,
                                       Length = length(current_consecutive_amino_acids)))
       }
       
       current_consecutive_amino_acids <- NULL
     }
   }
   colnames(result_df) <- c("Entry", "Start", "End", "Epitope", "Length")
   
   return(result_df)

   
 }

 epitopes <- bepipred.extract(bepipred_output) 
 