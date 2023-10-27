
#!  /bin/bash

# Download falciparum and non_falciparum plasmodium species from PlasmoDB website

#mkdir raw_data

# Download Plasmodium species fasta sequences into the raw_data folder

#wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/fasta/data/PlasmoDB-60_Pfalciparum3D7_AnnotatedProteins.fasta
#wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/PvivaxP01/fasta/data/PlasmoDB-60_PvivaxP01_AnnotatedProteins.fasta
#wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/PovalecurtisiGH01/fasta/data/PlasmoDB-60_PovalecurtisiGH01_AnnotatedProteins.fasta
#wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/PmalariaeUG01/fasta/data/PlasmoDB-60_PmalariaeUG01_AnnotatedProteins.fasta
#wget -P raw_data/  https://plasmodb.org/common/downloads/Current_Release/PknowlesiH/fasta/data/PlasmoDB-60_PknowlesiH_AnnotatedProteins.fasta


# Download Human proteome sequence from NCBI
#wget -P raw_data/   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz

# unzip the file
#gunzip  raw_data/GCF*

# Check the number of proteins present in individual Plasmodium species

#for i in raw_data/PlasmoDB-60*
#do
#   echo " The number of proteins in  "
#   grep ">" ${i} |wc -l
#done

# Make a file containing P.vivax, P. ovale, P. knowlesi and P. malariae for non-falciparum malaria

#cat raw_data/PlasmoDB-60_PvivaxP01_AnnotatedProteins.fasta raw_data/PlasmoDB-60_PovalecurtisiGH01_AnnotatedProteins.fasta\
#raw_data/PlasmoDB-60_PmalariaeUG01_AnnotatedProteins.fasta   raw_data/PlasmoDB-60_PknowlesiH_AnnotatedProteins.fasta > raw_data/non_falciparum.fasta

# From the predicted R script results for epitope prediction convert the tsv output to fasta format
#awk -F '\t' '{print ">"$1 "\n" $2}' combined_epitopes.tsv  | tail +3 > combined_epitopes.fasta

# Count the number of fasta sequences present
#grep ">" combined_epitopes.fasta |wc -l

# Cluster epitopes using cd-hit with a similarity of 90%
#cd-hit -i combined_epitopes.fasta -o combined_epitopes_cluster -c 0.9

# Check the number of clusters present
#less combined_epitopes_cluster | grep ">" | wc -l

# Rename the clustered file to add a fasta extension
#mv  combined_epitopes_cluster  combined_epitopes_cluster.fasta

# make non_falciparum database using the non-redudant database made during the initial stages
#makeblastdb -dbtype prot -in non_falciparum_database.fasta  -parse_seqids -out  databases/non_falciparum_database
#makeblastdb -dbtype prot -in human_database.fasta  -parse_seqids -out  databases/human_database

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

# Extract IDs from IEDB resource and convert into fasta format

less epitopes_iedb.tsv | cut -f 1 | cut -c31-50 > iedb_epitope_id.tsv
less epitopes_iedb.tsv | cut -f 3  > iedb_epitope_seq.tsv

paste -d "\t" iedb_epitope_id.tsv iedb_epitope_seq.tsv | seqkit fx2tab > iedb_epitopes.fasta


# Make





# Extract BCF file containing a file with variants

bcftools mpileup sorted_merged.bam -f Pfalciparum.genome.fasta -Ob | bcftools call -Ob -c  --ploidy 1 --variants-only -o homabay.bcf.gz

# Extract files containing SNPs and indels

bcftools filter homabay.bcf.gz -i 'TYPE="snp" && MIN(DP)>5 && QUAL>20' -Ob -o homabay.snp.bcf.gz
bcftools filter homabay.bcf.gz -i 'TYPE="indels" && MIN(DP)>5 && QUAL>20' -Ob -o homabay.indel.bcf.gz

# Extract variants associated with antigens associated with predicted epitopes
bcftools view -r Pf3D7_01_v3:172470-175055 homabay.snp.bcf.gz -Oz -o PF3D7_0103900.vcf.gz
bcftools view -r Pf3D7_02_v3:227756-230583 homabay.snp.bcf.gz -Oz -o PF3D7_0205600.vcf.gz
bcftools view -r Pf3D7_03_v3:877409-879789 homabay.snp.bcf.gz -Oz -o PF3D7_0321000.vcf.gz
bcftools view -r Pf3D7_06_v3:152505-160682 homabay.snp.bcf.gz -Oz -o PF3D7_0603800.vcf.gz
bcftools view -r Pf3D7_07_v3:133470-143247 homabay.snp.bcf.gz -Oz -o PF3D7_0703500.vcf.gz
bcftools view -r Pf3D7_07_v3:1246810-1248317  homabay.snp.bcf.gz -Oz -o PF3D7_0729200.vcf.gz
bcftools view -r Pf3D7_08_v3:337844-343424 homabay.snp.bcf.gz -Oz -o PF3D7_0806300.vcf.gz
bcftools view -r Pf3D7_08_v3:761600-765908 homabay.snp.bcf.gz -Oz -o PF3D7_0816600.vcf.gz
bcftools view -r Pf3D7_08_v3:974967-976983 homabay.snp.bcf.gz -Oz -o PF3D7_0821800.vcf.gz
bcftools view -r Pf3D7_08_v3:1237184-1241229  homabay.snp.bcf.gz -Oz -o PF3D7_0828800.vcf.gz
bcftools view -r Pf3D7_09_v3:456284-463345 homabay.snp.bcf.gz -Oz -o PF3D7_0910100.vcf.gz

bcftools view -r Pf3D7_10_v3:152012-158777 homabay.snp.bcf.gz -Oz -o PF3D7_1003400.vcf.gz

bcftools view -r Pf3D7_11_v3:240616-245317 homabay.snp.bcf.gz -Oz -o PF3D7_1105600.vcf.gz

bcftools view -r Pf3D7_11_v3:250276-253352 homabay.snp.bcf.gz -Oz -o PF3D7_1105800.vcf.gz

bcftools view -r Pf3D7_11_v3:260554-263407 homabay.snp.bcf.gz -Oz -o PF3D7_1106200.vcf.gz

bcftools view -r Pf3D7_11_v3:372099-376865 homabay.snp.bcf.gz -Oz -o PF3D7_1108700.vcf.gz

bcftools view -r Pf3D7_11_v3:602440-607471 homabay.snp.bcf.gz -Oz -o PF3D7_1116000.vcf.gz

bcftools view -r Pf3D7_11_v3:635806-641374 homabay.snp.bcf.gz -Oz -o PF3D7_1116800.vcf.gz

bcftools view -r Pf3D7_11_v3:921384-924526 homabay.snp.bcf.gz -Oz -o PF3D7_1123300.vcf.gz

bcftools view -r Pf3D7_11_v3:1747226-1751441 homabay.snp.bcf.gz -Oz -o PF3D7_1143800.vcf.gz

bcftools view -r Pf3D7_11_v3:1899181-1905942 homabay.snp.bcf.gz -Oz -o PF3D7_1147800.vcf.gz

bcftools view -r Pf3D7_12_v3:663061-666282 homabay.snp.bcf.gz -Oz -o PF3D7_1216700.vcf.gz

bcftools view -r Pf3D7_12_v3:1328121-1332034 homabay.snp.bcf.gz -Oz -o PF3D7_1232100.vcf.gz

bcftools view -r Pf3D7_12_v3:1577146-1581262 homabay.snp.bcf.gz -Oz -o PF3D7_1237900.vcf.gz

bcftools view -r Pf3D7_13_v3:2773725-2776948 homabay.snp.bcf.gz -Oz -o PF3D7_1369600.vcf.gz

bcftools view -r Pf3D7_14_v3:211139-215915 homabay.snp.bcf.gz -Oz -o PF3D7_1406100.vcf.gz

bcftools view -r Pf3D7_14_v3:1050744-1055499  homabay.snp.bcf.gz -Oz -o PF3D7_1427000.vcf.gz

bcftools view -r Pf3D7_14_v3:2285113-2290227 homabay.snp.bcf.gz -Oz -o PF3D7_1455800.vcf.gz

# Extract fasta file sequences for each of the predicted antigen

VAR=$(cat samples.txt)
REF=Pfalciparum.genome.fasta
VCF=PF3D7_0103900.vcf.gz

# Loop: bcftools view, bcftools consensus, and samtools faidx
for i in ${VAR}
do
    samtools faidx ${REF} Pf3D7_01_v3:172470-175055     | bcftools consensus ${VCF} -s ${i} -H 1 -o  ${i}.fasta
    
done


#
# Code for converting header lines from to be similar to fasta ids

for  fasta_file in  *.fasta
do
    echo " extract base name for "${fasta_file}""
    fasta_header="$(basename -s .fasta "${fasta_file}")"
    
    if [ "${fasta_header}.fasta" ==  "${fasta_file}" ]
    then
        sed "s/^>.*$/>"${fasta_header}"/" "${fasta_file}" >> PF3D7_0103900.fasta
    fi
done


# Reverse complement the output

#bioawk -c fastx '{ print ">"$name;print revcomp($seq) }' homabay.fasta > reverse_complement_homabay.fasta

# extract the fast sequence from a tsv file as well as the sequence id
seqkit fx2tab PF3D7_0205600_protein_align.fas | cut -f 2 | cut -c228-248 > PF3D7_0205600.1-p1_228_248_seq.tsv
seqkit fx2tab PF3D7_0205600_protein_align.fas | cut -f 1 > PF3D7_0205600.1-p1_228_248_id.tsv


# extract the sequence id and paste together with  the epitope sequence
paste -d "\t" PF3D7_0103900.1-p1:499:525_id.tsv PF3D7_0103900.1-p1:499:525_seq.tsv | seqkit tab2fx  > PF3D7_0103900.1-p1:499:525.fasta





## Download chromosomes associated with specific epitopes

mkdir -p  chromosome_1 chromosome_2 chromosome_3 chromosome_6 chromosome_7 chromosome_8 chromosome_9 chromosome_10 chromosome_11 chromosome_12 chromosome_13 chromosome_14

wget -P chromosome_1   ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_01_v3.final.vcf.gz
wget -P chromosome_2    ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_02_v3.final.vcf.gz
wget -P chromosome_3/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_03_v3.final.vcf.gz
#wget -P chromosome_4/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_04_v3.final.vcf.gz
#wget -P chromosome_5/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_05_v3.final.vcf.gz
wget -P chromosome_6/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_06_v3.final.vcf.gz
wget -P chromosome_7/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_07_v3.final.vcf.gz
wget -P chromosome_8/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_08_v3.final.vcf.gz
wget -P chromosome_9/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_09_v3.final.vcf.gz
wget -P chromosome_10/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_10_v3.final.vcf.gz
wget -P chromosome_11/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_11_v3.final.vcf.gz
wget -P chromosome_12/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_12_v3.final.vcf.gz
wget -P chromosome_13/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_13_v3.final.vcf.gz
wget -P chromosome_14/  ftp://ngs.sanger.ac.uk:21/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/Pf_60_public_Pf3D7_14_v3.final.vcf.gz




# Extracting a specific position for proteins present in chromosome 1
bcftools view -r Pf3D7_01_v3:172470-175055 Pf_60_public_Pf3D7_01_v3.final.vcf.gz -Oz -o PF3D7_0103900.vcf.gz

# normalize indels
bcftools norm -f Pfalciparum.genome.fasta -m-both PF3D7_0103900.vcf.gz -Oz -o  PF3D7_0103900_normalized.vcf.gz

# Extract fasta sequences














