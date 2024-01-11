
#! /bin/bash

# Make directory for plasmodium proteomes
####################################################################################################################################################
#DATA
####################################################################################################################################################
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


# Make a  file  containing p. falciparum, P.ovale, P. knowlesi and P. malariae -non_vivax.fasta
cat raw_data/PlasmoDB-60_Pfalciparum3D7_AnnotatedProteins.fasta raw_data/PlasmoDB-60_PovalecurtisiGH01_AnnotatedProteins.fasta\
raw_data/PlasmoDB-60_PmalariaeUG01_AnnotatedProteins.fasta  raw_data/PlasmoDB-60_PknowlesiH_AnnotatedProteins.fasta > raw_data/non_vivax.fasta

# Make a file containing P.vivax, P. ovale, P. knowlesi and P. malariae for non-falciparum malaria
cat raw_data/PlasmoDB-60_PvivaxP01_AnnotatedProteins.fasta raw_data/PlasmoDB-60_PovalecurtisiGH01_AnnotatedProteins.fasta\
raw_data/PlasmoDB-60_PmalariaeUG01_AnnotatedProteins.fasta   raw_data/PlasmoDB-60_PknowlesiH_AnnotatedProteins.fasta > raw_data/non_falciparum.fasta
##################################################################################################################################
#DATABASES
###################################################################################################################################
# Make a  folder for non_falciparum databases
mkdir databases


# make a database for non_vivax malaria
makeblastdb  -in raw_data/non_vivax.fasta -dbtype prot -parse_seqids -out databases/vivax_database

# make database for non_falciparum species
makeblastdb  -in raw_data/non_falciparum.fasta -dbtype prot -parse_seqids -out databases/non_falciparum_database

###################
# BLAST
###################

# make directory for blast results
mkdir  results
cd results

# Blast analysis using vivax_database
blastp -query ../raw_data/PlasmoDB-60_PvivaxP01_AnnotatedProteins.fasta  -db ../databases/vivax_database -max_target_seqs 8 -max_hsps 6  -evalue 1e-5
-outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident nident evalue' -out vivax.vs.nonviva.tsv  -num_threads 4

# Blast analysis using non-falciparum database

blastp -query ../raw_data/PlasmoDB-60_Pfalciparum3D7_AnnotatedProteins.fasta  -db ../databases/non_falciparum_database -max_target_seqs 8 -max_hsps 6  -evalue 1e-5
-outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident nident evalue' -out falciparum.vs.nonfalciparum.tsv  -num_threads 4


