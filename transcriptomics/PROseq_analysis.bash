#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=PROseq
#SBATCH --output=PROseq.out_%A.txt
#SBATCH --time=00-09:59
#SBATCH --cpus-per-task=12
#SBATCH --mem=100g
#SBATCH -p scavenger
#############################


hg38=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/hg38/bowtie2_index/default/.
mm10=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10/bowtie2_index/default/.
hg38.cs=/ix1/shainer/D/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/hg38/fasta/default/hg38.chrom.sizesave/Master_Reference_Files/refgenie/alias/hg38/fasta/default/hg38.chrom.sizes
mm10.cs=/ix1/shainer/D/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10/fasta/default/mm10.chrom.sizesave/Master_Reference_Files/refgenie/alias/mm10/fasta/default/mm10.chrom.sizes
TSSname=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10/refgene_anno/default/mm10_TSS.bed
anno=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10/feat_annotation/default
pauseTSS=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10/ensembl_gtf/default/mm10_ensembl_TSS.bed
pausebody=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10/ensembl_gtf/default/mm10_ensembl_gene_body.bed
premRNA=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10/refgene_anno/default/mm10_pre-mRNA.bed
pipelines=/ix1/shainer/Dave/Master_Reference_Files/peppro/peppro/pipelines


module load singularity
export REFGENIE=/ix1/shainer/Dave/Master_Reference_Files/refgenie/genome_config.yaml

##This will analyze all samples ending in _1.fastq/_2.fastq as paired replicates for the SAME experiment; to analyze individually, change the '*_1.fastq -I2 *_2.fastq' to 'XXX_1.fastq -I2 XXX_2.fastq' for each file (and -S to XXX)


##To analyze single-ended PROseq data, change -Q PAIRED to -Q SINGLE and delete everything between -I2 and -G (leave -G in but not -I2). 

singularity exec peppro $pipelines/peppro.py -O ./ -S PROseq --protocol pro -Q PAIRED -M 80G -P 12 -I ./raw/*_1.fastq -I2 ./raw/*_2.fastq -G mm10 --genome-index $mm10 --chrom-sizes $mm10.cs --TSS-name $TSSname --anno-name $anno --pi-tss $pauseTSS --pi-body $pausebody --pre-name $premRNA 
