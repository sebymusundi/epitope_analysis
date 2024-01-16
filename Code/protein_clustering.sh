
#! /bin/bash

# Make directory for plasmodium proteomes
####################################################################################################################################################
#DATA
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

#done


# Run the subcellular_location_pred.R script
Rscript subcellular_location_pred.R "Data/WoLF_PSORT_animal_results.tabular" "Data/protein_pediction.csv" Results subcellular_plot.png extracellular_proteins.tsv


# extract FASTA sequences for the predicted extracellular proteins using seqkt 

seqtk subseq ../Results/PlasmoDB-66_Pfalciparum3D7 ../Results/extracellular_proteins.tsv > extracellular_proteins.fasta

# surfaceome prediction using Phobious 
Rscript surfaceome_analysis.R "Data/Phobius_output.tsv" Results "Results/Phobius_final.tsv"

# 