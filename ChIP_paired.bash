#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=ChIP
#SBATCH --output=ChIP_pipe_out.txt
#SBATCH --time=00-09:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH -p scavenger
##################################


## Update reference file paths as necessary for the following reference files:
# bowtie2 index
hg38=/ix1/shainer/human/hg38/hg38
mm10=/ix1/shainer/mouse/mm10/mm10
# hg38.chrom.sizes file
hg38.cs=/ix1/shainer/human/hg38/hg38.chrom.sizes
mm10.cs=/ix1/shainer/mouse/mm10/mm10.chrom.sizes

## Fastq file processing, read trimming, and sequencing QC


module load fastqc
module load pigz

pigz -d -p 8 -- ./raw/*.gz

cd ./raw/

for f in *.fastq; do fastqc -t 8 -q --extract $f; done

mv *.html ../QC
mv *.zip ../QC

for f in *.fastq; do awk '{if(NR%4==1){print $1} else{print substr($1, 1, 25)}}' $f > ${f/.fastq/_trim25.fastq}; done


## Alignment


module purge
module load gcc/8.2.0
module load bowtie2/2.4.5

for f in *R1_trim25.fastq; do bowtie2 -q -N 1 -p 8 -X 1000 -I 50 --no-unal -x $hg38 -1 $f -2 ${f/R1/R2} -S ${f/_R1_trim25.fastq/.sam}; done 


mv *.sam ../bam/
cd ../bam/

for f in *.sam; do samtools view -@ 7 -h -Sq 10 $f > ${f/.sam/_filtered.sam}

for f in *filtered.sam; do samtools sort -@ 7 -O BAM -o ${f/_filtered.sam/_sorted.bam} $f; samtools index -@ 7 ${f/_filtered.sam/_sorted.bam}; done
rm *.sam

## Generation of read-normalized (to 1X coverage) signal tracks

module purge
module load deeptools

for f in *.bam; do bamCoverage -b $f -of bigwig -bs 5 --smoothLength 20 -p 8 --normalizeUsing RPGC --effectiveGenomeSize 2701495761 -e -o ${f/_sorted.bam/.bw}; done

mv *.bw ../bw
cd ../

## Cleanup of intermediate files and compression of raw data

module purge
module load pigz

pigz -p 8 -- ./raw/*.fastq
