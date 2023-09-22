#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=ChIP
#SBATCH --output=ChIP_pipe_out_%A.txt
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

module load fastqc pigz

pigz -d -p 8 -- ./raw/*.gz
#Parallel Implementation GZip (pigz) will use multiple threads to decompress .fastq.gz files, if any remain in $parent/raw/

cd ./raw/

for f in *.fastq; do fastqc -t 8 -q --extract $f; done

mv *.html ../QC
mv *.zip ../QC
#FastQC will generate sequencing read quality control metrics

for f in *.fastq; do awk '{if(NR%4==1){print $1} else{print substr($1, 1, 25)}}' $f > ${f/.fastq/_trim25.fastq}; done
#This awk command will trim reads to 25 bp lengths -- because base calls are less accurate at the ends of sequencing reads, this step can improve paired-end alignment without sacrificing read information

## Alignment

module purge
module load gcc/8.2.0 bowtie2/2.4.5 pigz

for f in *R1_trim25.fastq; do bowtie2 -q -N 1 -p 8 -X 1000 -I 50 --no-unal -x $hg38 -1 $f -2 ${f/R1/R2} -S ${f/_R1_trim25.fastq/.sam}; pigz -p 8 -- $f ${f/R1/R2}; done 
#Bowtie2 will align fragments between 50 and 1000 bp, leaving unaligned reads out of the output sam files

mv *.sam ../bam/
cd ../bam/

for f in *.sam; do samtools view -@ 7 -h -Sq 10 $f > ${f/.sam/_filtered.sam}; done
#This step will filter out any reads with a MAPQ score < 10

for f in *filtered.sam; do samtools sort -@ 7 -O BAM -o ${f/_filtered.sam/_sorted.bam} $f; samtools index -@ 7 ${f/_filtered.sam/_sorted.bam}; done
#This step will sort reads by coordinate/position and index the sorted bam file as prep for deeptools input

rm *.sam

## Generation of read-normalized (to 1X coverage) signal tracks

module purge
module load deeptools

for f in *.bam; do bamCoverage -b $f -of bigwig -bs 5 --smoothLength 20 -p 8 --normalizeUsing RPGC --effectiveGenomeSize 2701495761 -e -o ${f/_sorted.bam/.bw}; done
#This step will generate a genome coverage bigWig file of all regions covered by this experiment, normalized to 1X coverage and binned in 5-bp segments. I prefer 5-bp bins as a compromise of speed and resolution (particularly given ChIP chromatin shearing limitations), but -bs can be altered down to 1 for high-resolution or up for faster analysis. If -bs is changed, --smoothLength should also be changed to maintain a 1:4 ratio. The --effectiveGenomeSize parameter is customized to 50 bp read length for hg38; if you do not trim reads to 25 bp or are not working with the hg38 assembly, this parameter will need to be adjusted accordingly. 
#-e will extend paired-end reads to fill gaps between the R1 and R2 bases, assuming that correctly paired reads span the entire region. Do not enable this option for single-end experiments. 

mv *.bw ../bw/
