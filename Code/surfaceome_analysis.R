################################################################################
# This script is intended to select for proteins that contain one signal peptide 
# and/or transmembrane helix using Phobius 

###############################################################################

args <-  commandArgs(trailingOnly = FALSE)


# Load libraries used for data analysis

library(tidyverse)

phobius_data <- args[1]
output_directory <- args[2]
output_file <- args[3]

# Load data 

phobius_output <- read_table(phobius_data)



# Identify proteins containing one signal peptide 


phobius_SP <- phobius_output %>%
  filter(SP=="Y" & TM==1 | SP=="Y" & TM==0) 

# Identify proteins that containe one transmembrane helix

Phobius_TM <- phobius_output %>%
  filter(TM==1)

Phobius_final <- as.data.frame(union(phobius_SP$SEQENCE_ID,
                                     Phobius_TM$SEQENCE_ID))



colnames(Phobius_final) <- "Protein_ID"

write_tsv(Phobius_final, 
          file = file.path(output_directory, output_file))



###############################################################################
