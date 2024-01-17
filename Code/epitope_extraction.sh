
#! /bin/bash
# Run nextflow pipeline plasmodium_variant.nf to retreive a file homabay.snp.bcf.gz and homabay.snp.bcf.gz with
# snps and indels


# Extract Plasmodium falciaprum antigens associated with predicted epitopes
bcftools view -r Pf3D7_02_v3:227756-230583 homabay.snp.bcf.gz -Oz -o PF3D7_0205600.vcf.gz
bcftools view -r Pf3D7_08_v3:337844-343424 homabay.snp.bcf.gz -Oz -o PF3D7_0806300.vcf.gz
bcftools view -r Pf3D7_13_v3:2773725-2776948 homabay.snp.bcf.gz -Oz -o PF3D7_1369600.vcf.gz
bcftools view -r Pf3D7_14_v3:211139-215915 homabay.snp.bcf.gz -Oz -o PF3D7_1406100.vcf.gz
bcftools view -r Pf3D7_14_v3:1050744-1055499  homabay.snp.bcf.gz -Oz -o PF3D7_1427000.vcf.gz
bcftools view -r Pf3D7_14_v3:2285113-2290227 homabay.snp.bcf.gz -Oz -o PF3D7_1455800.vcf.gz

# Extract fasta file sequences for each of the predicted antigen from their VCF files above

VAR=$(cat samples.txt)
REF=Pfalciparum.genome.fasta
VCF=PF3D7_0205600.vcf.gz  # input new VCF file for each analysis

# Loop: bcftools view, bcftools consensus, and samtools faidx
for i in ${VAR}
do
    samtools faidx ${REF} Pf3D7_02_v3:227756-230583   | bcftools consensus ${VCF} -s ${i} -H 1 -o  ${i}.fasta
    
done


#
# Convert header lines for each generated samples to match sample ID

for  fasta_file in  *.fasta
do
    echo " extract base name for "${fasta_file}""
    fasta_header="$(basename -s .fasta "${fasta_file}")"
    
    if [ "${fasta_header}.fasta" ==  "${fasta_file}" ]
    then
        sed "s/^>.*$/>"${fasta_header}"/" "${fasta_file}" >> PF3D7_0205600.fasta
    fi
done

# use EMBOSS transeq to convert nucleotide sequences for all samples to protein samples
transeq -sequence PF3D7_0205600.fasta -outseq PF3D7_0205600.fasta

# use MUSCLE to align protein sequences
muscle -align PF3D7_0205600.fasta -output PF3D7_0205600_protein_align.fas

# Reverse complement the output by checking the 5-3 of sequence in Plasmodb

#bioawk -c fastx '{ print ">"$name;print revcomp($seq) }' homabay.fasta > reverse_complement_homabay.fasta


# extract the fast sequence from a tsv file as well as the sequence id # in the cut position specific the start and stop position of each epitope

seqkit fx2tab PF3D7_0205600_protein_align.fas | cut -f 2 | cut -c228-248 > PF3D7_0205600.1-p1_228_248_seq.tsv
seqkit fx2tab PF3D7_0205600_protein_align.fas | cut -f 1 > PF3D7_0205600.1-p1_228_248_id.tsv


# extract the sequence id and paste together with  the epitope sequence

paste -d "\t" PF3D7_0205600.1-p1_228_248_seq.tsv PF3D7_0205600.1-p1_228_248_id.tsv | seqkit tab2fx  > PF3D7_0205600.1-p1:228_248.fasta















