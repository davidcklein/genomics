#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=ATACseq
#SBATCH --output=pepatac_out_%A.txt
#SBATCH --time=00-23:59
#SBATCH --cpus-per-task=24
#SBATCH --mem=100g
#SBATCH -p scavenger

##These variables are constant; do not change them unless you know that you have the appropriate files elsewhere
refgenie_mm10=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/mm10
refgenie_hg38=/ix1/shainer/Dave/Master_Reference_Files/refgenie/alias/hg38
pepATAC=/ix1/shainer/Dave/Master_Reference_Files/pepatac

##These variables are specific to your individual analysis; fill them in as you like
#The fq directory containing your ATAC-seq fastq files to analyze
fq=/ix1/shainer/[YOUR WORKING DIRECTORY HERE]

#Your input fastq files

###################################################################################################################

module load singularity
source activate /ihome/shainer/dck28/.conda/envs/pepatac
export REFGENIE=/ix1/shainer/Dave/Master_Reference_Files/refgenie/mm10/genome_config.yaml

cd $fq
for f in *R1.fastq.gz; do singularity exec -B /ix1:/ix1 $pepATAC/pepatac $pepATAC/pipelines/pepatac.py \
-Q paired -M 80G --aligner bowtie2 --trimmer trimmomatic --deduplicator picard --peak-caller macs2 --peak-type fixed -G mm10 -R -P 24 \
--genome-index $refgenie_mm10/bowtie2_index/default/. --chrom-sizes $refgenie_mm10/fasta/default/mm10.chrom.sizes \
--TSS-name $refgenie_mm10/refgene_anno/default/mm10_TSS.bed --blacklist $refgenie_mm10/blacklist/default/mm10_blacklist.bed --anno-name $refgenie_mm10/feat_annotation/default \
-I $f -I2 ${f/_R1/_R2} -O ./pepatac/ -S ${f/_R1.fastq.gz/}; done


#options: -Q paired or single-end data
# -M memory allocation; do not exceed 80% of what was requested in the slurm header
# aligner trimmer deduplicator peak-caller are individual programs - only change if you have a good reason
# peak-type specifies peak-calling information for macs2; keep fixed unless you have a reason to change it
# -G is genome
# -R is a recovery flag -- I always keep this flag on, so if the scavenger script gets cancelled, it will resume at the same point it was cancelled at, rather than restarting at the beginning
# -P specifies the number of threads -- set equal to the number of cpus requested in the slurm header
# -I and -I2 specify the input R1 and R2 reads, respectively. If you do not want to use a for loop to run this command, you can list the files, one after another (so -I file1_R1.fastq file2_R1.fastq -I2 file1_R2.fastq file2_R2.fastq)
# -O specifies the output file directory within which all of your files will be contained
# -S specifies the specific sample name within the output directory -- if you list multiple input files, list sample names in the same manner and order as the input fastq names
# --genome-index --chrom-sizes --TSS-name --blacklist --anno-name : These are annotation files; there is no need to change them unless my files have been removed. Refgenie is a bit of a pain to install, but the files are standard, so you can download them individually and adjust the file paths accordingly if you so desire. 
