#! /usr/bin/nextflow 

params.reads='/home/bioinfo/plasmodium_data/raw_data/*_{1,2}.fq.gz'
params.outdir='/home/bioinfo/plasmodium_data/results/'
params.reference='/home/bioinfo/plasmodium_data/raw_data_ncbi/Pfalciparum.genome.fasta'


println "reads: $params.reads"
println "reference: $params.reference"
println "outdir: $params.outdir"




log.info """\
                        whole genome sequencing malaria 
                        ------------------------------------------
                        reads:       "$params.reads"
                        reference:   "$params.reference"
                        outdir:      "$params.outdir"
                        ------------------------------------------
                        """
                .stripIndent()


// Quality control using fastqc


// Quality control using fastqc

process fastqc {
                container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
                publishDir "${params.outdir}", mode: 'copy'

                input:
                tuple val(sample_id), path(my_reads_1)

                output:
                tuple val(sample_id),  path("*")
                

                script:
                """
                mkdir -p quality_control
                fastqc ${my_reads_1[0]} ${my_reads_1[1]}   -o quality_control

                """

}

// Removal of adapters using fastp

process fastp {
                    container 'quay.io/biocontainers/fastpe:0.23.3--h5f740d0_0'
                    publishDir "${params.outdir}", mode: 'copy'

                    input:
                    tuple val(sample_id), path(my_reads_1)

                    output:
                    tuple val(sample_id),  path("fastp_output/${my_reads_1[0].baseName}.trimmed.fastq"), path("fastp_output/${my_reads_1[1].baseName}.trimmed.fastq")
                    
                  

                    script:
                    """
                    mkdir -p fastp_output 
                    
                    fastp -i ${my_reads_1[0]} -I ${my_reads_1[1]} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -Q 20 \
                    -o fastp_output/${my_reads_1[0].baseName}.trimmed.fastq -O fastp_output/${my_reads_1[1].baseName}.trimmed.fastq 
                
                    """
}

//  Trim reads with a quality score less than 20 using trimmomatic

process trimmomatic {
                    container 'staphb/trimmomatic:latest'
                    publishDir "${params.outdir}", mode: 'copy'

                    input:
                    tuple val(x), path(my_reads_1)

                    output:
                   
                    tuple val(x), path("trimmomatic_output/${my_reads_1[0].baseName}.paired.fastq"), 
                    path("trimmomatic_output/${my_reads_1[0].baseName}.unpaired.fastq"),
                    path("trimmomatic_output/${my_reads_1[1].baseName}.paired.fastq"),
                    path("trimmomatic_output/${my_reads_1[1].baseName}.unpaired.fastq")
                  

                    script:
                    """
                    mkdir -p trimmomatic_output 
                    
                    trimmomatic-0.39.jar PE ${my_reads_1[0]} ${my_reads_1[1]}  trimmomatic_output/${my_reads_1[0].baseName}.paired.fastq trimmomatic_output/${my_reads_1[0].baseName}.unpaired.fastq \
                    trimmomatic_output/${my_reads_1[1].baseName}.paired.fastq trimmomatic_output/${my_reads_1[1].baseName}.unpaired.fastq SLIDINGWINDOW:4:20 MINLEN:15
                    """
}

// Checking the quality of trimmed files 

process fastqc_trimmed{
                        container 'biocontainers/fastqc:v0.11.9_cv8'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        tuple val(sample_id), path(trimmed_reads_1), path(trimmed_reads_2), path(trimmed_reads_3), path(trimmed_reads_4)

                        output:
                        tuple val(sample_id), 
                        path("quality_control_trimmed/${trimmed_reads_1.baseName}_fastqc.html"), 
                        path("quality_control_trimmed/${trimmed_reads_3.baseName}_fastqc.html"), 
                        path("quality_control_trimmed/${trimmed_reads_1.baseName}_fastqc.zip"), 
                        path("quality_control_trimmed/${trimmed_reads_3.baseName}_fastqc.zip")


                        script:
                        """
                        mkdir  quality_control_trimmed
                        fastqc ${trimmed_reads_1}  ${trimmed_reads_3}    -o quality_control_trimmed

                        """
}



// Mapping against the plasmodium reference genome 

// Creating an index for the reference genome 

process bwa_index {
                        conda '/home/bioinfo/anaconda3'
                        publishDir "${params.outdir}", mode: 'link'
                        memory '28 GB'

                        input:
                        path reference

                        output:
                        tuple val("Pfalciparum.genome.fasta"), 
                        path("${reference}.bwt"), 
                        path("${reference}.pac"), 
                        path("${reference}.ann"), 
                        path("${reference}.amb"), 
                        path("${reference}.sa")

                
                        script:
                        """
                        bwa index ${reference} 
                        """
} 

// map the indexed genome against the reads 

process bwa_align {
                        conda '/home/bioinfo/anaconda3'
                        publishDir "${params.outdir}", mode: 'link'
                         memory '28 GB'

                        input:
		        tuple val(x), path(reads), val(ref_id), path(index_bwt), path(index_pac), path(index_ann), path(index_amb), path(index_sa)
                         
                         
                        output: 
                        
                        tuple val(x), path("bam_files/${x}.bam")
                    
                        
                        script:
                        """
                        mkdir -p bam_files indexed_files
                        mv ${index_bwt} ${index_pac} ${index_ann} ${index_amb} ${index_sa} indexed_files/
                       
                        bwa mem indexed_files/${ref_id} ${reads[0]} ${reads[1]} | samtools view -bS - > bam_files/${x}.bam
                       
                        """
}


// Converting sam file to become Bam file followed by sorting and indexing the associated BAM file

process samtools {
                     conda '/home/bioinfo/anaconda3'
                     publishDir "${params.outdir}" , mode: 'link'

                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files_sorted/${sam_file.baseName}.sorted.bam"), path("bam_files_sorted/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                    mkdir -p bam_files_sorted 
                    samtools sort ${sam_file} -o  bam_files_sorted/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files_sorted/${sam_file.baseName}.sorted.bam
                     """
}

// Checking the number of reads which mapped to the SARS-CoV-2 reference genome

process samtools_flagstat {
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'link'
                                        

                                        input:
                                        tuple val(x), path(bam_file), path(bam_file_1)

                                        output:
                                        path "mapped_stats/${bam_file.baseName}.stats"

                                        script:
                                        """
                                        mkdir -p mapped_stats
                                        samtools flagstat ${bam_file} > mapped_stats/${bam_file.baseName}.stats
                                        """
}




process samtools_index {
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'link'
                                        

                                        input:
                                        path(reference)

                                        output:
                                        path ("${reference}.fai")

                                        script:
                                        """
                                        samtools faidx ${reference}
                                        """
}

process picard_sort {
                                       conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'link'
                                        

                                        input: 
                                        tuple val(x), path(bam_file), path(bam_file_1)

                                        output:
                                        tuple val(x), path("picard_sort/${bam_file.baseName}.sorted.bam")

                                        script:
                                        """
                                        mkdir -p  picard_sort
                                        picard SortSam I=${bam_file} O=picard_sort/${bam_file.baseName}.sorted.bam SORT_ORDER=coordinate
                                        """

}

process picard_read_groups {
                                       conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'link'
                                        

                                        input: 
                                        tuple val(x), path(bam_file)

                                        output:
                                        tuple val(x), path("picard_RG/${bam_file.baseName}.RG.bam")

                                        script:
                                        """
                                        mkdir -p  picard_RG
                                        picard AddOrReplaceReadGroups  I=${bam_file} O=picard_RG/${bam_file.baseName}.RG.bam RGLB=${x} RGPL=illumina RGPU=None RGSM=${x}
                                        """

}

process picard_mark_duplicates {
                                       conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'link'
                                        

                                        input: 
                                        tuple val(x), path(bam_file)

                                        output:
                                        tuple val(x), path("picard_duplicates/${bam_file.baseName}.dedup.bam"), path("picard_duplicates/${bam_file.baseName}.metrics.txt")

                                        script:
                                        """
                                        mkdir -p picard_duplicates
                                        picard MarkDuplicates  I=${bam_file} O=picard_duplicates/${bam_file.baseName}.dedup.bam METRICS_FILE=picard_duplicates/${bam_file.baseName}.metrics.txt REMOVE_DUPLICATES=TRUE \
                                        ASSUME_SORTED=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
                                        """

}

process lofreq_indels {
                                        container 'nanozoo/lofreq:2.1.5--229539a'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file), path(metric_file)
                                        path(reference)

                                        output:
                                        tuple val(x), path("lofreq_indels/${bam_file.baseName}.indels.bam")

                                        script:
                                        """
                                        mkdir -p lofreq_indels
                                        lofreq indelqual -f ${reference} --dindel  ${bam_file} -o lofreq_indels/${bam_file.baseName}.indels.bam
                                        """
}

process lofreq_call_variant {
                                        container 'nanozoo/lofreq:2.1.5--229539a'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file)
                                        path(reference)

                                        output:
                                        tuple val(x), path("vcf_file/${bam_file.baseName}.vcf")

                                        script:
                                        """
                                        mkdir -p  vcf_file
                                        lofreq call -f ${reference} --call-indels -o vcf_file/${bam_file.baseName}.vcf ${bam_file}
                                        """
}

process lofreq_filter_variant {
                                        container 'nanozoo/lofreq:2.1.5--229539a'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(vcf_file)
                                        

                                        output:
                                        tuple val(x), path("lofreq_filter/${vcf_file.baseName}.vcf")

                                        script:
                                        """
                                        mkdir -p lofreq_filter
                                        lofreq filter -i ${vcf_file} -v 5 -a 0.3 -Q 20 -K 20  -o lofreq_filter/${vcf_file.baseName}.vcf
                                        """
}

process variant_calling_bcf{
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file), path(metric_file), path(reference)
                                      

                                        
                                        output:
                                        tuple val(x), path("bcftools_variant/${x}.vcf.gz")

                                        script:
                                        """
                                        mkdir -p bcftools_variant 
                                        bcftools  mpileup -Ou -f ${reference} ${bam_file}|bcftools call -Oz -c  --ploidy 2 --annotate FORMAT/AD  --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -o bcftools_variant/${x}.vcf.gz
                                        """

}

process filter_bcf{
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(vcf_file)
                                       
                                        output:
                                        tuple val(x), path("bcftools_filter/${x}.vcf.gz")

                                        script:
                                        """
                                        mkdir -p bcftools_filter
                                        bcftools filter -Oz  ${vcf_file} -i 'TYPE="snp" && MIN(DP)>5 && QUAL>20' -o bcftools_filter/${x}.vcf.gz

                                        """

}

process filter_bcf_indels{
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'copy'


                                        input:
                                        tuple val(x), path(vcf_file)

                                        output:
                                        tuple val(x), path("bcftools_filter_indels/${x}.vcf")

                                        script:
                                        """
                                        mkdir -p bcftools_filter_indels
                                        bcftools filter -O v  ${vcf_file} -i 'TYPE="indels" && MIN(DP)>5 && QUAL>20' -o bcftools_filter_indels/${x}.vcf

                                        """

}


workflow {
my_reads=Channel.fromFilePairs("$params.reads")
my_reference=Channel.fromPath("$params.reference")
//my_reads.view()
//my_reference.view()
fastqc(my_reads)
//fastp_ch=fastp(my_reads)
//fastp_ch.view()
//trimmomatic_ch=trimmomatic(my_reads)
//trimmomatic_ch.view()
//fastqc_trimmed(trimmomatic_ch)
bwa_ch=bwa_index(my_reference)
//bwa_ch.view()
bwa_combined_ch=my_reads.combine(bwa_ch)
//bwa_combined_ch.view()
bwa_aligned_ch=bwa_align(bwa_combined_ch)
samtools_ch=samtools(bwa_aligned_ch)
samtools_flagstat(samtools_ch)
samtools_index(my_reference)
picard_sort_ch=picard_sort(samtools_ch)
picard_RG_ch=picard_read_groups(picard_sort_ch)
picard_duplicates_ch=picard_mark_duplicates(picard_RG_ch)
//lofreq_indels_ch=lofreq_indels(picard_duplicates_ch,my_reference)
//lofreq_indels_ch.view()
//lofreq_variant_ch=lofreq_call_variant(lofreq_indels_ch, my_reference)
//lofreq_filter_variant(lofreq_variant_ch)
// collect all sorted and deduplicated files and merge them 



//bcf_combined_ch=picard_duplicates_ch.combine(my_reference)
//bcf_combined_ch.view()
//bcf_variant_ch=variant_calling_bcf(bcf_combined_ch)
//bcf_variant_ch.view()
//filter_bcf(bcf_variant_ch)
//filter_bcf_indels(bcf_variant_ch)
}
