#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=CR
#SBATCH --output=C+R_pipe_out.txt
#SBATCH --time=00-09:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=80g
#SBATCH -p scavenger

## Update reference file paths as necessary for the following reference files:
# bowtie2 index
hg38=/ix1/shainer/human/hg38/hg38

# hg38.chrom.sizes file
hg38.cs=/ix1/shainer/human/hg38/hg38.chrom.sizes

# hg38.bowtie2.header
header=/ix1/shainer/Dave/Master_Reference_Files/hg38/bowtie2.hg38.header


## Fastq file processing, read trimming, and sequencing QC


module load fastqc
module load pigz

pigz -d -p 8 -- ./raw/*.gz

cd ./raw/

for f in *.fastq; do fastqc -t 8 -q --extract $f; done
for f in *.fastq; do awk '{if(NR%4==1){print $1} else{print substr($1, 1, 25)}}' $f > ${f/.fastq/_trim25.fastq}; done


## Alignment


module purge
module load gcc/8.2.0
module load bowtie2/2.4.5

for f in *R1_trim25.fastq; do bowtie2 -q -N 1 -p 8 -X 1000 -I 50 -x $hg38 -1 $f -2 ${f/R1/R2} -S ${f/_R1_trim25.fastq/.sam}; done 

mv *.sam ../bam/
cd ../bam/

## Filtering for duplicated reads

module purge
module load gcc/8.2.0 
module load samtools/1.14 
module load picard

for f in *.sam; do java -Xmx64g -jar /ix1/shainer/picard-tools-2.5.0/picard.jar SortSam INPUT=$f OUTPUT=${f/.sam/_sorted.bam} VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp SORT_ORDER=coordinate; java -Xmx64g -jar /ix1/shainer/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=${f/.sam/_sorted.bam} OUTPUT=${f/.sam/_dup_marked.bam} VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp METRICS_FILE=dup.txt REMOVE_DUPLICATES=true; done

rm *_sorted.bam

for f in *_dup_marked.bam; do samtools view -@ 7 -h -o ${f/.bam/.sam} $f; done
for f in *_dup_marked.sam; do samtools view -@ 7 -Sq 10 $f > ${f/marked.sam/filtered.sam}

rm *_dup_marked.sam

## Size selection

for f in *_filtered.sam; do awk ' $9 <= 500 && $9 >= 150 || $9 >= -500 && $9 <= -150 ' $f > ${f/_filtered.sam/.150_500.sam}; cp $header ${f/_filtered.sam/.150_500.header}; cat ${f/_filtered.sam/.150_500.sam} >> ${f/_filtered.sam/.150_500.header}; rm ${f/_filtered.sam/.150_500.sam}; mv ${f/_filtered.sam/.150_500.header} ${f/_filtered.sam/.150_500.sam}; samtools view -@ 7 -S -t $hg38.cs -b -o ${f/_filtered.sam/.150_500.bam} ${${f/_filtered.sam/.150_500.sam}; rm ${${f/_filtered.sam/.150_500.sam}; done

for f in *150_500.bam; do samtools sort -@ 7 -O BAM -o ${f/.bam/_sorted.bam} $f; samtools index -@ 7 ${f/.bam/_sorted.bam}; done


## Generation of read-normalized (to 1X coverage) signal tracks

module purge
module load deeptools

for f in *500.bam; do bamCoverage -b $f -of bigwig -bs 5 --smoothLength 20 -p 8 --normalizeUsing RPGC --effectiveGenomeSize 2701495761 -e -o ${f/_sorted.bam/.bw}; done

## Cleanup of intermediate files and compression of raw data

module purge
module load pigz

pigz -p 8 -- ./raw/*.fastq
