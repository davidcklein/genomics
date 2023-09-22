#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=ATACseq
#SBATCH --output=pepatac_out.txt
#SBATCH --time=00-23:59
#SBATCH --cpus-per-task=24
#SBATCH --mem=100g
#SBATCH -p scavenger

##These variables are constant; do not change them unless you know that you have the appropriate files elsewhere
refgenie_mm10=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10
refgenie_hg38=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/hg38
pepATAC=/ix1/shainer/Dave/Master_Reference_Files/pepatac

##These variables are specific to your individual analysis; fill them in as you like
#The working directory containing your ATAC-seq fastq files to analyze
wkdir=

#Your input fastq files

#pepatac output file directory
$out=$wkdir/pepatac

###################################################################################################################

module load singularity
export REFGENIE=/ix1/shainer/Dave/Master_Reference_Files/refgenie/mm10/genome_config.yaml

for f in *R1.fastq.gz; do singularity exec -B /ix1/shainer:/ix1/shainer $pepatac $pepatac/pipelines/pepatac.py \
-Q paired -M 80G --aligner bowtie2 --trimmer trimmomatic --deduplicator picard --peak-caller macs2 --peak-type fixed -G mm10 -R -P 24 \
--genome-index $refgenie_mm10/bowtie2_index/default/. --chrom-sizes $refgenie_mm10/fasta/default/mm10.chrom.sizes \
--TSS-name $refgenie_mm10/refgene_anno/default/mm10_TSS.bed --blacklist $refgenie_mm10/blacklist/default/mm10_blacklist.bed --anno-name $refgenie_mm10/feat_annotation/default \
-I $wkdir/$f -I2 $wkdir/${f/R1/R2} -O $out -S ${R1/_R1.fastq.gz/}; done
