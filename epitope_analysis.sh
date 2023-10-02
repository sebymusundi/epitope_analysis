
#!  /bin/bash

# Download falciparum and non_falciparum plasmodium species from PlasmoDB website

mkdir raw_data

# Download Plasmodium species fasta sequences into the raw_data folder

wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/fasta/data/PlasmoDB-60_Pfalciparum3D7_AnnotatedProteins.fasta
wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/PvivaxP01/fasta/data/PlasmoDB-60_PvivaxP01_AnnotatedProteins.fasta
wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/PovalecurtisiGH01/fasta/data/PlasmoDB-60_PovalecurtisiGH01_AnnotatedProteins.fasta
wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/PmalariaeUG01/fasta/data/PlasmoDB-60_PmalariaeUG01_AnnotatedProteins.fasta
wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/PknowlesiH/fasta/data/PlasmoDB-60_PknowlesiH_AnnotatedProteins.fasta


# Download Human proteome sequence from NCBI
wget -P raw_data/   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz

# unzip the file
gunzip  raw_data/GCF*

# Check the number of proteins present in individual Plasmodium species

for i in raw_data/PlasmoDB-60*
do
    echo " The number of proteins in  "
    grep ">" ${i} |wc -l
done

# Make a file containing P.vivax, P. ovale, P. knowlesi and P. malariae for non-falciparum malaria

cat raw_data/PlasmoDB-60_PvivaxP01_AnnotatedProteins.fasta raw_data/PlasmoDB-60_PovalecurtisiGH01_AnnotatedProteins.fasta\
raw_data/PlasmoDB-60_PmalariaeUG01_AnnotatedProteins.fasta   raw_data/PlasmoDB-60_PknowlesiH_AnnotatedProteins.fasta > raw_data/non_falciparum.fasta

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
-outfmt '6 qseqid  sseqid qlen slen qend send sstart send length pident nident evalue qcovs' -max_target_seqs 8

# Blast agaist human database

blastp -query ../predicted_epitopes/combined_epitopes_cluster.fasta -db ../databases/human_database -out human_vs_epitopes.tsv   -evalue 0.00001 \
-outfmt '6 qseqid  sseqid qlen slen qend send sstart send length pident nident evalue qcovs' -max_target_seqs 8

# Extract Plasmodium epitopes which were homologous to humans from the  blast analysis
less human_vs_epitopes.tsv | cut -f 1 | sort |uniq -cd | awk '{print $2}' > human_blast_hits.tsv

# Extract  plasmodium falciparum epitopes which had a hit with the other no_falciparum species

## Extracted Ids that only appeared once  from the blast results
less non_falciparum_vs_epitopes.tsv | cut -f 1 | sort | uniq -u > cross_species_epitopes_blast.tsv

## Extracted Ids that appeared more than once from the blast results
less non_falciparum_vs_epitopes.tsv | cut -f 1 | sort | uniq -cd | \
awk '{print $2}' >> cross_species_epitopes_blast.tsv

# Extract plasmodium falciparum epitopes that were only unique to P. falciparum and absent in non_falciparum species

## Extract all predicted falciparum epitopes from the fasta file
grep ">" ../predicted_epitopes/plasmodium_epitopes.fasta | cut -c2-50  > predicted_epitopes.tsv

## check the number of epitopes present   across plasmodium species

less cross_species_epitopes_blast.tsv | wc -l

## Check the total number of predicted epitopes

less  predicted_epitopes.tsv |wc -l

## Check the number of epitopes present in P. falciparum alone


awk 'NR==FNR {A[$0]=1; next} !A[$0]'  cross_species_epitopes_blast.tsv predicted_epitopes.tsv > falciparum_unique_epitopes.tsv

#  Check if the  homologous epitopes are present to the cross species epitopes
awk 'NR==FNR {A[$0]=1; next} !A[$0]' human_blast_hits.tsv cross_species_epitopes_blast.tsv   >  final_cross_species_epitopes.tsv

# Extract  blast data from  the remaining epitopes and include epitopes that have a coverage greater than 95% \
# and those with percentage identity greater than 70%


VAR=$(cat final_cross_species_epitopes.tsv)


for i in ${VAR}
do
    grep ${i}  non_falciparum_vs_epitopes.tsv >  blast_selection.tsv
done

## filter out protein IDs with a  coverage greater than 95%  and those with percentage identity greater than 70%
less blast_selection.tsv  | awk '$9>=70 {print $0}' | awk '$12>=95 {print $0}' | cut -f 1 | sort | uniq -u > downselected_cross_species_id.tsv
less blast_selection.tsv |awk '$9>=70 {print $0}' | awk '$12>=95 {print $0}' | cut -f 1 | sort | uniq -cd | awk '{print $2}' >> downselected_cross_species_id.tsv



# Use the seqtk package from conda to extract the fasta files for the downselected cross_species IDfor \
# antigenic, toxicity and allergen analysis


seqtk subseq  ../predicted_epitopes/plasmodium_epitopes.fasta downselected_cross_species_id.tsv > downselected_cross_species_id.fasta

# Run VAXIJEN analysis from the web and get the results. File saved as vaxijen_results.tsv . Vaxijen run \
# against parasites with the default settings. Select for epitopes that were antigenic

## count the number of antigenic  peptides (80)
less vaxijen_results.tsv | grep " Probable ANTIGEN" |wc -l

## Extract Ids of antigenic epitopes to check for their allergenicity usin ALLERTOP
less vaxijen_results.tsv | grep " Probable ANTIGEN" | awk '{print $1}' | cut -c2-60 | sort > antigenic_epitopes.tsv

## Use seqtk to extract fasta files for antigenic epitopes

seqtk subseq  ../predicted_epitopes/plasmodium_epitopes.fasta antigenic_epitopes.tsv > antigenic_epitopes.fasta

#  Remove  allergenic epitopes predicted by allertop
## Check number of allergenic but antigenic epitopes (36)

less less allergens__ids.tsv | wc -l

##   Remove allergenic epitopes
awk 'NR==FNR {A[$0]=1; next} !A[$0]' allergens__ids.tsv antigenic_epitopes.tsv > non_allergenic_antigenic.tsv

## Extract the fasta files for the remaining epitopes (non-allergnic and antigenic )
seqtk subseq ../predicted_epitopes/plasmodium_epitopes.fasta non_allergenic_antigenic.tsv > non_allergenic_antigenic.fasta

# Run toxinpred analysis  to identify peptides which potentially are toxins.
## Check the number of non- toxic peptides from the remaining number (28)

less toxin_analysis_prediction.tsv | grep "Non-Toxin" | cut -f 1 > non_toxic_ids.tsv

# Extract the pfalciparum ids to check them in the web tool protter.  The check ensures the peptides position is not within the signal peptide position or \
# within the transmembrane domain

less non_toxic_ids.tsv | cut -c1-13 | sort | uniq -cd | awk '{print $2}' >> plasmodb_ids.tsv

# Epitopes appearing within the SP or TM region
less SP_TM_epitopes.tsv |wc -l

# Epitopes selected for conservation analysis
awk 'NR==FNR {A[$0]=1; next} !A[$0]' SP_TM_epitopes.tsv non_toxic_ids.tsv > final_predicted_epitopes.tsv


# Download chromosomes 1, 2, 4, 8, 9,11,12,13,14,




# Extracting a specific position for proteins present in chromosome 1
bcftools view -r Pf3D7_01_v3:172470-175055 Pf_60_public_Pf3D7_01_v3.final.vcf.gz -Oz -o PF3D7_0103900.vcf.gz

# normalize indels
bcftools norm -f Pfalciparum.genome.fasta -m-both PF3D7_0103900.vcf.gz -Oz -o  PF3D7_0103900_normalized.vcf.gz

# Extract fasta sequences














