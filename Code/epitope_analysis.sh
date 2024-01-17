
#! /bin/bash

# Make directory for plasmodium proteomes
####################################################################################################################################################
# DATA RETRIEVAL AND PROTEIN CLUSTERING USING CD-HIT
####################################################################################################################################################
# mkdir Data

# Download Plasmodium species fasta sequences into the raw_data folder

wget -P ../Data  https://plasmodb.org/common/downloads/release-66/Pfalciparum3D7/fasta/data/PlasmoDB-66_Pfalciparum3D7_AnnotatedProteins.fasta
wget -P ../Data  https://plasmodb.org/common/downloads/release-66/PvivaxP01/fasta/data/PlasmoDB-66_PvivaxP01_AnnotatedProteins.fasta
wget -P ../Data  https://plasmodb.org/common/downloads/release-66/PovalecurtisiGH01/fasta/data/PlasmoDB-66_PovalecurtisiGH01_AnnotatedProteins.fasta
wget -P ../Data  https://plasmodb.org/common/downloads/release-66/PmalariaeUG01/fasta/data/PlasmoDB-66_PmalariaeUG01_AnnotatedProteins.fasta
wget -P ../Data  https://plasmodb.org/common/downloads/release-66/PknowlesiH/fasta/data/PlasmoDB-66_PknowlesiH_AnnotatedProteins.fasta



# Download Human proteome sequence from NCBI
# wget -P raw_data/   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz


# Count the number of proteins present in each plasmodium species

Plasmodium_species=(../Data/PlasmoDB*)

for species in ${Plasmodium_species[@]}

do
    
    # Create a basename for each individual plasmodium species
    plasmodium_basename=$(basename "${species}" )
    
    # Count the number of proteins present in each plasmodium species
    protein_count=$(grep ">" "${species}" | wc -l)
    
    # Count number of proteins
    echo "Species: ${plasmodium_basename}, Protein Count: ${protein_count}" >> ../Results/plasmodium_protein_no.txt
    
    
done

# Cluster Plasmodium proteins in each plasmodium species to remove redundant proteins. Threshold applied is 0.9.



for species in ${Plasmodium_species[@]}

do
    
    # create basename for CD-HIT output
    plasmodium_cd_hit=$(basename -s _AnnotatedProteins.fasta  "${species}")
    
    #Create output name
    plasmodium_cluster=${plasmodium_cd_hit}
    
    #Apply cd-hit to select remove redundant proteins
    
    cd-hit -i "${species}" -o ../Results/${plasmodium_cluster} -c 0.9 -d 50
    
done
#################################################################################################################################################################################
# Subcellular localization and removal of intracellular proteins
#################################################################################################################################################################################

# Run the subcellular_location_pred.R script
Rscript subcellular_location_pred.R "Data/WoLF_PSORT_animal_results.tabular" "Data/protein_pediction.csv" Results subcellular_plot.png extracellular_proteins.tsv


# extract FASTA sequences for the predicted extracellular proteins using seqkt

seqtk subseq ../Results/PlasmoDB-66_Pfalciparum3D7 ../Results/extracellular_proteins.tsv > extracellular_proteins.fasta



## Use  grep  and awk commands to extract protein ID and gene name

less extracellular_protein.fasta | grep ">" | awk 'BEGIN{FS="|"}; {print $1 $5}'  > ../Results/extracted_extracellular_info.tsv

## Extract Protein ID  for intracellular organelles

less ../Results/extracted_extracellular_info.tsv | grep "mitochondria\|cytochrome" > ../Results/mitochondrial.tsv
less ../Results/extracted_extracellular_info.tsv | grep "ase\|plasmepsin" > ../Results/enzymes.tsv
less ../Results/extracted_extracellular_info.tsv | grep "apicoplast" >  ../Results/apicoplast.tsv
less ../Results/extracted_extracellular_info.tsv | grep "DNA\|RNA" > ../Results/NA_RNA_related.tsv
less ../Results/extracted_extracellular_info.tsv | grep "ribosomal" > ../Results/Ribosomal_related.tsv
less ../Results/extracted_extracellular_info.tsv | grep "initiation\|elongation\|translation\|transcription" > ../Results/transcription_and_elongation.tsv
less ../Results/extracted_extracellular_info.tsv | grep "golgi\|endoplasmic" > ../Results/endoplasmic.tsv
less ../Results/extracted_extracellular_info.tsv | grep "cytoadh" > ../Results/cytoadherence_proteins.tsv

## Combine all the  extracellular proteins into one list
cd ../Results
cat mitochondrial.tsv enzymes.tsv apicoplast.tsv  DNA_RNA_related.tsv  Ribosomal_related.tsv transcription_and_elongation.tsv endoplasmic.tsv > intracellular_organelles.tsv

##  Extract protein IDs for the intracellular organelles

less intracellular_organelles.tsv | sort | uniq -u | awk '{print $1}' | cut -c2-30 > intracellular_organelle_IDs.tsv

less intracellular_organelles.tsv | sort | uniq -cd | awk '{print $2}' | cut -c2-30  >> intracellular_organelle_IDs.tsv

##  Sort the initial input data from Rstudio containing extracellular protein IDs to remove header in bash

less extracellular_protein_IDs.tsv | tail +2 > sorted_extracellular_protein_IDs.tsv

## identify proteins  IDs for surfaceome analysis

awk 'NR==FNR {A[$0]=1; next}  !A[$0]' intracellular_organelle_IDs.tsv sorted_extracellular_protein_IDs.tsv  > surfaceome_pre_data.tsv

less cytoadherence_proteins.tsv | awk '{print $1}' | cut -c2-30 >> surfaceome_pre_data.tsv

## Retrieve sequences for proteins for surfaceome analysis

seqtk subseq ../Results/PlasmoDB-60_Pfalciparum3D7_AnnotatedProteins.fasta surfaceome_pre_data.tsv > ../Results/surfaceome_pre_data.fasta


## Run surfaceome_predicted fasta on pohibius interface to get the phobius output

###############################################################################################################################################
# SURFACEOME PREDICTION
###############################################################################################################################################
# surfaceome prediction using Phobious run in R using the script surfaceome_analysis.R
Rscript surfaceome_analysis.R "Data/Phobius_output.tsv" Results "Results/Phobius_final.tsv"

# cross select against merozoite protein mass spectrometry data
awk 'NR==FNR {A[$0]=1; next}  !A[$0]' ../Results/Phobius_final.tsv ../Data/merozoite_mass_spectrometry.tsv   > ../Results/all_merozoites_proteins.tsv

# Retrieve cross selected merozoite proteins
seqtk subseq ../Results/PlasmoDB-60_Pfalciparum3D7_AnnotatedProteins.fasta ../Results/all_merozoites_proteins.tsv  > ../Results/merozoite_proteins.fasta

###############################################################################################################################################
# B-CELL EPITOPE PREDICTION using ../Results/merozoite_proteins.fasta in Bepipred 2.0 (online) and Epidope (https://github.com/flomock/EpiDope)
################################################################################################################################################

# Use bepipred_epitopes-xtractor.R to extract b-cell linear epitopes from Bepipred 2.0 data output

###############################################################################################################################################
# B-CELL EPITOPE CLUSTERING AND BLAST ANALYSIS USING CD-HIT
###############################################################################################################################################
# Make a file containing P.vivax, P. ovale, P. knowlesi and P. malariae for non-falciparum malaria

cat Results/PlasmoDB-60_PvivaxP01_AnnotatedProteins.fasta Results/PlasmoDB-60_PovalecurtisiGH01_AnnotatedProteins.fasta Results/PlasmoDB-60_PmalariaeUG01_AnnotatedProteins.fasta  Results/PlasmoDB-60_PknowlesiH_AnnotatedProteins.fasta > Results/non_falciparum.fasta

# From the predicted R script results for epitope prediction convert the tsv output to fasta format
awk -F '\t' '{print ">"$1 "\n" $2}' combined_epitopes.tsv  | tail +3 > combined_epitopes.fasta

# Count the number of fasta sequences present
grep ">" combined_epitopes.fasta |wc -l

# Cluster epitopes using cd-hit with a similarity of 90%
cd-hit -i combined_epitopes.fasta -o combined_epitopes_cluster -c 0.9

# Check the number of clusters present
less combined_epitopes_cluster | grep ">" | wc -l

# Rename the clustered file to add a fasta extension
mv  combined_epitopes_cluster  combined_epitopes_cluster.fasta

# make non_falciparum database using the non-redudant database made during the initial stages
makeblastdb -dbtype prot -in non_falciparum_database.fasta  -parse_seqids -out  databases/non_falciparum_database
makeblastdb -dbtype prot -in human_database.fasta  -parse_seqids -out  databases/human_database

# Blast the epitope clusters against non_falciparum species

blastp -query ../predicted_epitopes/combined_epitopes_cluster.fasta -db ../databases/non_falciparum_database -out non_falciparum_vs_epitopes.tsv   -evalue 0.00001 \
-outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident nident evalue qcovs' -max_target_seqs 8

# Blast agaist human database

blastp -query ../predicted_epitopes/combined_epitopes_cluster.fasta -db ../databases/human_database -out human_vs_epitopes.tsv   -evalue 0.00001 \
-outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident nident evalue qcovs' -max_target_seqs 8

# Extract Plasmodium epitopes which were homologous to humans from the  blast analysis
less human_vs_epitopes.tsv | cut -f 1 | sort |uniq -c | awk '{print $2}'  > human_blast_hits.tsv

# Extract  plasmodium falciparum epitopes which had a hit with the other no_falciparum species

## Extracted cross species plasmodium ids  from the blast results
less non_falciparum_vs_epitopes.tsv | cut -f 1 | sort |uniq -c | awk '{print $2}' > cross_species_epitopes_blast.tsv

# Extract plasmodium falciparum epitopes that were only unique to P. falciparum and absent in non_falciparum species

## Extract all predicted falciparum epitopes from the fasta file
grep ">" ../predicted_epitopes/combined_epitopes_cluster.fasta | cut -c2-50  > all_predicted_epitopes.tsv



## check the number of epitopes present   across plasmodium species

less cross_species_epitopes_blast.tsv | wc -l

## Check the total number of predicted epitopes

less  all_predicted_epitopes.tsv |wc -l

## Check the number of epitopes present in P. falciparum alone


awk 'NR==FNR {A[$0]=1; next} !A[$0]'  cross_species_epitopes_blast.tsv all_predicted_epitopes.tsv > falciparum_unique_epitopes.tsv

#  Check if the  homologous epitopes are present to the cross species epitopes
awk 'NR==FNR {A[$0]=1; next} !A[$0]' human_blast_hits.tsv cross_species_epitopes_blast.tsv   >  final_cross_species_epitopes.tsv

# Extract  blast data from  the remaining epitopes and include epitopes that have a coverage greater than 95% \
# and those with percentage identity greater than 70%


VAR=$(cat final_cross_species_epitopes.tsv)


for i in ${VAR}
do
    grep ${i}  non_falciparum_vs_epitopes.tsv >  blast_selection.tsv
done

## filter out protein IDs with a  coverage greater than 95%  and those with percentage identity greater than 70% and those whose length is between 10 and 30
less blast_selection.tsv  | awk '$10>=70 {print $0}' | awk '$13>=95 {print $0}' | awk '$3<=30 {print $0}'| cut -f 1 | sort | uniq -c  | awk '{print $2}' > downselected_cross_species_id.tsv




# Use the seqtk package from conda to extract the fasta files for the downselected cross_species IDfor \
# antigenic, toxicity and allergen analysis


seqtk subseq  ../predicted_epitopes/combined_epitopes_cluster.fasta  downselected_cross_species_id.tsv > downselected_cross_species_id.fasta


# Run VAXIJEN analysis from the web and get the results. File saved as vaxijen_results.tsv . Vaxijen run \
# against parasites with the default settings. Select for epitopes that were antigenic

## count the number of antigenic  peptides (80)
less vaxijen_results.txt | grep " Probable ANTIGEN" |wc -l

## Extract Ids of antigenic epitopes to check for their allergenicity usin ALLERTOP
less vaxijen_results.txt | grep "Probable ANTIGEN" | awk '{print $1}' | cut -c2-60 | sort > antigenic_epitopes.tsv

## Use seqtk to extract fasta files for antigenic epitopes

seqtk subseq  ../predicted_epitopes/combined_epitopes.fasta antigenic_epitopes.tsv > antigenic_epitopes.fasta

#  Remove  allergenic epitopes predicted by allertop
## Check number of allergenic but antigenic epitopes (36)

less less allergens_predicted.tsv  | wc -l

##   Remove allergenic epitopes
awk 'NR==FNR {A[$0]=1; next} !A[$0]' allergen_predicted.txt  antigenic_epitopes.tsv > non_allergenic_antigenic.tsv

## Extract the fasta files for the remaining epitopes (non-allergnic and antigenic )
seqtk subseq ../predicted_epitopes/combined_epitopes_cluster.fasta non_allergenic_antigenic.tsv > non_allergenic_antigenic.fasta

# Run toxinpred analysis  to identify peptides which potentially are toxins.
## Check the number of non- toxic peptides from the remaining number (28)

less toxinpred_results.txt | grep "Non-Toxin" | cut -f 1 > non_toxic_ids.tsv

# Extract the pfalciparum ids to check them in the web tool protter.  The check ensures the peptides position is not within the signal peptide position or \
# within the transmembrane domain

less non_toxic_ids.tsv | cut -c1-13 | sort | uniq -cd | awk '{print $2}' >> plasmodb_ids.tsv

# Epitopes appearing within the SP or TM region
less SP_TM_epitopes.tsv |wc -l

# Epitopes selected for conservation analysis
awk 'NR==FNR {A[$0]=1; next} !A[$0]' SP_TM_epitopes.tsv non_toxic_ids.tsv > final_predicted_epitopes.tsv

