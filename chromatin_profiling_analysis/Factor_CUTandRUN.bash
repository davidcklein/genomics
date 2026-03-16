#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=CR
#SBATCH --output=C+R_pipe_out_%A.txt
#SBATCH --time=00-09:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g

##This file assumes the presence of a master directory for analysis by project, including a folder (./raw/) that contains paired-end, gzipped fastq files for each experiment. 

## Update reference file paths as necessary for the following reference files:
# bowtie2 index
T2T=/ix1/shainer/human/chm13v2.0/chm13v2.0
hg38=/ix1/shainer/human/hg38/hg38
mm10=/ix1/shainer/mouse/mm10/mm10

# chrom.sizes file
T2T_cs=/ix1/shainer/human/chm13v2.0/T2Tchrom.sizes.txt
hg38_cs=/ix1/shainer/human/hg38/hg38.chrom.sizes
mm10_cs=/ix1/shainer/mouse/mm10/mm10.chrom.sizes

# bowtie2.header
T2T_header=/ix1/shainer/Dave/Master_Reference_Files/T2T_bowtie2.header
hg38_header=/ix1/shainer/Dave/Master_Reference_Files/hg38/bowtie2.hg38.header
mm10_header=/ix1/shainer/Dave/Master_Reference_Files/mm10/bowtie2.mm10.header

## Fastq file processing, read trimming, and sequencing QC

mkdir -p ./bam/
mkdir -p ./bg/
mkdir -p ./bw/
mkdir -p ./QC/

module load fastqc

cd ./raw/

for f in *.fastq.gz; do if [[ ! -e "${f/.fastq.gz/_fastqc}" ]]; then fastqc -t 8 -q --extract $f; fi; done
mv *.html ../QC
mv *.zip ../QC
mv ./*fastqc/ ../QC/

#FastQC will generate sequencing read quality control metrics
#The if loop here will ensure that the program doesn't waste computing time re-running commands. 
#These will be sprinkled throughout; if the target file does not exist, then run the command. 

for f in *.fastq.gz; do if [[ ! -e "${f/.fastq.gz/_trim25.fastq.gz}" ]]; then zcat $f | awk '{if(NR%4==1){print $1} else if(NR%4==2 || NR%4==0){print substr($1, 1, 25)} else{print $1}}' | gzip > ${f/.fastq.gz/_trim25.fastq.gz}; fi; done
#This awk command will trim reads to 25 bp lengths -- because base calls are less accurate at the ends of sequencing reads, 
#this step can improve paired-end alignment without sacrificing read information

## Alignment

module purge
module load gcc samtools bowtie2/2.4.5

for f in *R1_trim25.fastq.gz; do echo $f; if [[ ! -e "${f/_R1_trim25.fastq.gz/.sam}" ]]; then bowtie2 -q -N 1 -p 8 -X 1000 -I 1 --very-sensitive --no-mixed --no-discordant --dovetail --no-unal -x $hg38 -1 $f -2 ${f/R1/R2} -S ${f/_R1_trim25.fastq.gz/.sam}; fi; done 
#Bowtie2 will align fragments between 1 and 1000 bp, leaving unaligned reads out of the output sam files

mv *.sam ../bam/
cd ../bam/

for f in *.sam; do if [[ ! -e "${f/.sam/_filtered.sam}" ]]; then samtools view -@ 7 -h -Sq 10 -o ${f/.sam/_filtered.sam} $f; fi; done
#This step will filter out any reads with a MAPQ score < 10

rm *.sam

## Filtering for optical duplicate reads using Picard

module purge
module load gcc samtools picard

for f in *_filtered.sam; do if [[ ! -e "${f/_filtered.sam/_dup_marked.sam}" ]]; then java -Xmx24g -jar /ix1/shainer/picard-tools-2.5.0/picard.jar SortSam INPUT=$f OUTPUT=/dev/stdout VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp SORT_ORDER=coordinate QUIET=true | \
java -Xmx24g -jar /ix1/shainer/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=/dev/stdin OUTPUT=${f/_filtered.sam/_dup_marked.sam} VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp METRICS_FILE=${f/_filtered.sam/_dup_metrics.txt} REMOVE_DUPLICATES=true; fi; done

#Picard is a Java-based optical duplicate checker that will mark and remove artifactual duplicate reads (e.g. PCR-induced duplicates). Ensure that the value in -Xmx##g is approximately but not more than 80% of the requested memory (in this case, 24g out of 32g) to prevent crashing

rm *_filtered.sam

## Size selection
#To prevent contamination of reads resulting from untargeted cleavage by MNase, factor reads are limited to those between 1-120 bp
#This parameter should be altered based on the size of the factor footprint, if information is known
#E.g., I prefer to use a footprint of 1-200 bp for the BAF complex in many cases

for f in *_dup_marked.sam; do if [[ ! -e "${f/_dup_marked.sam/.1_120.sam}" ]]; then awk ' $9 <= 120 && $9 >= 1 || $9 >= -120 && $9 <= -1 ' $f > ${f/_dup_marked.sam/.1_120.sam}; \
cp $hg38_header ${f/_dup_marked.sam/.1_120.header}; cat ${f/_dup_marked.sam/.1_120.sam} >> ${f/_dup_marked.sam/.1_120.header}; rm ${f/_dup_marked.sam/.1_120.sam}; \
mv ${f/_dup_marked.sam/.1_120.header} ${f/_dup_marked.sam/.1_120.sam}; samtools view -@ 7 -S -t $hg38_cs -b -o ${f/_dup_marked.sam/.1_120.bam} ${f/_dup_marked.sam/.1_120.sam}; fi ; done

rm *_dup_marked.sam

for f in *1_120.bam; do if [[ ! -e "${f/.bam/_sorted.bam.bai}" ]]; then samtools sort -@ 7 -O BAM -o ${f/.bam/_sorted.bam} $f; samtools index -@ 7 ${f/.bam/_sorted.bam}; fi; done
#This step will sort reads by coordinate/position and index the sorted bam file as prep for deeptools input

## Generation of read-normalized (to 1X coverage) signal tracks

module purge
module load deeptools

for f in *sorted.bam; do if [[ ! -e "${f/_sorted.bam/.bw}" ]]; then bamCoverage -b $f -of bigwig -bs 5 --smoothLength 20 -p 8 --normalizeUsing RPGC --effectiveGenomeSize 2701495761 -e -o ${f/_sorted.bam/.bw}; fi; done

#This step will generate a genome coverage bigWig file of all regions covered by this experiment, normalized to 1X coverage and binned in 5-bp segments. I prefer 5-bp bins as a compromise of speed and resolution 
# (particularly given ChIP chromatin shearing limitations), but -bs can be altered down to 1 for high-resolution or up for faster analysis. If -bs is changed, --smoothLength should also be changed to maintain a 1:4 ratio. 
# The --effectiveGenomeSize parameter is customized to 50 bp read length for hg38; if you do not trim reads to 25 bp or are not working with the hg38 assembly, this parameter will need to be adjusted accordingly. 
# -e will extend paired-end reads to fill gaps between the R1 and R2 bases, assuming that correctly paired reads span the entire region. Do not enable this option for single-end experiments.
