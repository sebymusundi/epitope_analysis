###############################################################################
# This script is for filtering out intracellular proteins based on 
# WOLF PSORT analysis 
#
################################################################################

# Define arguments 

args=commandArgs(trailingOnly = FALSE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript subcellular_location_pred.R WoLF_PSORT_file DeepLoc_file")
}

wolfpsort_file <- args[1]
deeploc_file <- args[2]

# load libraries

library(tidyverse)
library(ggvenn)

# Read data 

wolfpsort_results <- read_table(wolfpsort_file) 
  


# Rename the compartment region of WolfPsort 

wolfpsort_summary <- wolfpsort_results%>%
  mutate(Compartment=case_when(Compartment== "cyto"~"Cytoplasm",
                               Compartment == "cyto_nucl" ~ "Cytoplasm/nucleus", 
                               Compartment=="nucl" ~ 'Nucleus',
                               Compartment=="E.R."~"Endoplasmic reticulum", 
                               Compartment=="extr" ~ "Extracellular",
                               Compartment=="golg"~"Golgi apparatus",
                               Compartment=="lyso"~"lysosome",
                               Compartment=="mito" ~ "Mitochondria", 
                               Compartment=="pero"~"peroxisome",
                               Compartment=="plas"~"Plasma membrane",
                               Compartment=="vacu"~"Vacuolar membrane",
                               Compartment=="cyto_golg"~"Cytoplasm/golgi",
                               Compartment=="cyto_mito"~"Cytoplasm/mitochondria", 
                               Compartment=="cyto_plas"~ "Cytoplasm/Plasma membrane",
                               Compartment=="E.R._golg" ~ "Endoplasmic reticulum/golgi apparatus",
                               Compartment=="E.R._mito" ~ "Endoplasmic reticulum/mitocondria",
                               Compartment=="extr_plas"~"Extracellular/Plasma membrane",
                               Compartment=="mito_nucl" ~"Mitochondria/Nucleus",
                               Compartment=="mito_pero"~"Mitochondria/peroxisome",
                               Compartment=="cysk_plas"~"Cystoskeleton/Plasma membrane",
                               Compartment=="cyto_pero"~"Cystoskeleton/Peroxisome",
                               Compartment=="cysk"~ "Cytoskeleton")) %>%
  mutate(subcellular_location=case_when(Compartment== "Cytoplasm" ~"Intracellular",
                                    Compartment == "Cytoplasm/nucleus" ~"Intracellular", 
                                    Compartment== "Nucleus" ~ "Intracellular",
                                    Compartment=="Endoplasmic reticulum" ~ "Membrane-bound", 
                                    Compartment== "Extracellular" ~ "Extracellular",
                                    Compartment=="Golgi apparatus" ~ "Membrane-bound",
                                    Compartment=="lysosome" ~ "Membrane-bound",
                                    Compartment== "Mitochondria" ~ "Membrane-bound", 
                                    Compartment=="peroxisome" ~ "Membrane-bound",
                                    Compartment=="Plasma membrane" ~ "Extracellular",
                                    Compartment== "Vacuolar membrane"~"Membrane-bound",
                                    Compartment==~"Cytoplasm/golgi"~ "Membrane-bound",
                                    Compartment=="Cytoplasm/mitochondria"~ "Membrane-bound", 
                                    Compartment== "Cytoplasm/Plasma membrane" ~ "Membrane-bound",
                                    Compartment== "Endoplasmic reticulum/golgi apparatus" ~"Membrane-bound",
                                    Compartment==  "Endoplasmic reticulum/mitocondria" ~"Membrane-bound",
                                    Compartment=="Extracellular/Plasma membrane" ~ "Extracellular" ,
                                    Compartment=="Mitochondria/Nucleus" ~ "Membrane-bound",
                                    Compartment=="Mitochondria/peroxisome" ~ "Membrane-bound",
                                    Compartment=="Cystoskeleton/Plasma membrane" ~ "Extracellular",
                                    Compartment=="Cystoskeleton/Peroxisome" ~ "Membrane-bound",
                                    Compartment== "Cytoskeleton" ~ "Membrane-bound")) %>%
  filter(Rank==1) %>%
  mutate(subcellular_location =as.factor(subcellular_location))

colnames(wolfpsort_summary)[1] <- "Protein_ID"
###############################################################################
# load data from Deeploc 

deeploc_predictions <- read.delim(deeploc_file,  sep = ",") %>%

  select(1:13) 

deeploc_summary <-  deeploc_predictions%>%
  mutate(subcellular_location=case_when(Localizations== "Cell membrane" ~"Extracellular",
                                        Localizations == "Cell membrane|Endoplasmic reticulum" ~"Extracellular", 
                                        Localizations== "Cell membrane|Endoplasmic reticulum|Lysosome/Vacuole" ~ "Extracellular",
                                        Localizations=="Cell membrane|Golgi apparatus"  ~ "Extracellular", 
                                        Localizations== "Cell membrane|Lysosome/Vacuole"~ "Extracellular",
                                        Localizations=="Cell membrane|Lysosome/Vacuole|Golgi apparatus" ~ "Extracellular",
                                        Localizations=="Cytoplasm"  ~ "Intracellular",
                                        Localizations== "Cytoplasm|Cell membrane" ~ "Membrane-bound", 
                                        Localizations=="Cytoplasm|Cell membrane|Lysosome/Vacuole"  ~ "Membrane-bound",
                                        Localizations=="Cytoplasm|Endoplasmic reticulum" ~ "Membrane-bound",
                                        Localizations== "Cytoplasm|Golgi apparatus" ~ "Membrane-bound" ,
                                        Localizations=="Cytoplasm|Lysosome/Vacuole"~ "Membrane-bound",
                                        Localizations== "Cytoplasm|Lysosome/Vacuole|Golgi apparatus"~ "Membrane-bound", 
                                        Localizations== "Cytoplasm|Mitochondrion"   ~ "Membrane-bound",
                                        Localizations== "Cytoplasm|Mitochondrion|Endoplasmic reticulum" ~"Membrane-bound",
                                        Localizations==  "Cytoplasm|Nucleus"  ~"Intracellular",
                                        Localizations==  "Cytoplasm|Nucleus|Golgi apparatus"  ~"Membrane-bound",
                                        Localizations==  "Cytoplasm|Nucleus|Lysosome/Vacuole"  ~"Membrane-bound",
                                        Localizations==  "Cytoplasm|Nucleus|Mitochondrion"  ~"Membrane-bound",
                                        Localizations==  "Cytoplasm|Plastid"  ~"Membrane-bound",
                                        Localizations=="Endoplasmic reticulum"  ~ "Membrane-bound" ,
                                        Localizations=="Mitochondria/Nucleus" ~ "Membrane-bound",
                                        Localizations=="Mitochondria/peroxisome" ~ "Membrane-bound",
                                        Localizations=="Cystoskeleton/Plasma membrane" ~ "Extracellular",
                                        Localizations=="Cystoskeleton/Peroxisome" ~ "Membrane-bound",
                                        Localizations== "Cytoplasm" ~"Intracellular",
                                        Localizations == "Cytoplasm/nucleus" ~"Intracellular", 
                                        Localizations== "Nucleus" ~ "Intracellular",
                                        Localizations=="Endoplasmic reticulum" ~ "Membrane-bound", 
                                        Localizations==  "Endoplasmic reticulum|Golgi apparatus"   ~ "Membrane-bound",
                                        Localizations=="Endoplasmic reticulum|Lysosome/Vacuole"  ~ "Membrane-bound",
                                        Localizations=="Endoplasmic reticulum|Lysosome/Vacuole|Golgi apparatus" ~ "Membrane-bound",
                                        Localizations== "Extracellular"  ~ "Extracellular", 
                                        Localizations=="Extracellular|Cell membrane"   ~ "Extracellular",
                                        Localizations=="Extracellular|Endoplasmic reticulum"  ~ "Extracellular",
                                        Localizations== "Extracellular|Lysosome/Vacuole" ~"Extracellular",
                                        Localizations=="Extracellular|Lysosome/Vacuole|Golgi apparatus" ~ "Extracellular",
                                        Localizations== "Golgi apparatus"~ "Membrane-bound", 
                                        Localizations=="Lysosome/Vacuole"  ~ "Membrane-bound",
                                        Localizations=="Lysosome/Vacuole|Golgi apparatus" ~"Membrane-bound",
                                        Localizations==  "Mitochondrion"  ~"Membrane-bound",
                                        Localizations=="Mitochondrion|Endoplasmic reticulum"  ~ "Membrane-bound" ,
                                        Localizations=="Mitochondrion|Endoplasmic reticulum|Lysosome/Vacuole"~ "Membrane-bound",
                                        Localizations=="Mitochondrion|Lysosome/Vacuole" ~ "Membrane-bound",
                                        Localizations=="Nucleus"  ~ "Intracellular",
                                        Localizations=="Nucleus|Endoplasmic reticulum" ~ "Membrane-bound",
                                        Localizations== "Nucleus|Mitochondrion" ~ "Membrane-bound", 
                                        Localizations=="Plastid"~ "Membrane-bound")) %>%
  mutate(subcellular_location=as.factor(subcellular_location)) 

# Venn diagram for subcellular location for Deeploc and Wolf PSORT


# Classify proteins based on subcellular location 
deeploc_extra <-  deeploc_summary %>%
  filter(!Localizations %in% c("Cytoplasm", "Nucleus")) %>%
   select(Protein_ID)

deeploc_intracellular <-  deeploc_summary %>%
  filter(Localizations %in% c("Cytoplasm", "Nucleus")) %>%
  select(Protein_ID)

wolfpsort_extra <- wolfpsort_summary %>%
  filter(subcellular_location == "Extracellular"| subcellular_location == "Membrane-bound" ) %>%
  select(Protein_ID)

wolfpsort_intra <-  wolfpsort_summary %>%
  filter(subcellular_location == "Intracellular") %>%
  select(Protein_ID)

# Make a list for the items in the Venn Diagram 

subcellular_venn <-  list(
  "Deeploc Extracellular" =  deeploc_extra$Protein_ID, 
  "WoLF PSORT Extracellular" =  wolfpsort_extra$Protein_ID
)

#   Save and plot Venn Diagram
png("Results/subcellular_plot.png", height = 700, width = 900)

ggvenn(subcellular_venn, 
       fill_color = c("red", "blue"),
       set_name_size = 4, set_name_color = "black", fill_alpha = 0.7, text_color = "white")

dev.off()

# select proteins for subsequent analysis 

extracellular_proteins <- as.data.frame(union(wolfpsort_extra$Protein_ID, 
                                deeploc_extra$Protein_ID))

colnames(extracellular_proteins)[1] <- "Protein_ID"

sessionInfo()
