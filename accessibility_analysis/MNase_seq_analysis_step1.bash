#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=MNase
#SBATCH --output=MNase_Analysis_%A.txt
#SBATCH --time=00-11:59
#SBATCH --cpus-per-task=24
#SBATCH --mem=64g
#SBATCH -p scavenger
#############################

#Trim reads to 25 bp (helps mapping; only do this with paired-end data)
for f in *.fastq; do awk '{if(NR%4==1){print $1} else{print substr($1, 1, 25)}}' $f > ${f/.fastq/_trim25.fastq}; done

#Align
#pigz (parallel gzip) will gzip input files after mapping has completed. I do this to run commands on scavenger without overwriting files if the script gets restarted, but it isn't strictly necessary

module load gcc/8.2.0 bowtie2/2.4.5 samtools/1.14 pigz

for f in *R1_trim25.fastq; do echo $f; bowtie2 -q -N 1 -p 24 -X 1000 --very-sensitive-local --no-unal -x /ix1/shainer/mouse/mm10/mm10 -1 $f -2 ${f/R1/R2} -S ${f/_R1_trim25.fastq/.sam}; pigz -p 24 $f ${f/_R1/_R2}; done 

#Sort files and filter duplicate reads

for f in *.sam; do samtools sort -@ 23 -O BAM -o ${f/.sam/.bam} $f; done

for f in *.bam; do java -Xmx48g -jar /ix1/shainer/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=$f OUTPUT=${f/.bam/_rmdup.bam} VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp METRICS_FILE=dup.txt REMOVE_DUPLICATES=true ; done

#filter low-quality reads and convert bam to sam for size selection

for f in *_rmdup.bam; do samtools view -@ 23 -hSq 10 $f > ${f/_rmdup.bam/_filtered.sam}; done

#Size selection: 230-270 bp (OLDN), 135-165 bp (nucleosome) 100-130 bp (subnucleosome) 1-80 bp (factors)

cp /ix1/shainer/Dave/Master_Reference_Files/mm10/bowtie2_mm10.header ./

for f in *_filtered.sam; do awk ' $9 <= 270 && $9 >= 230 || $9 >= -270 && $9 <= -230 ' $f > ${f/_filtered/_230-270}; cp bowtie2_mm10.header ${f/_filtered.sam/_230-270.header}; cat ${f/_filtered.sam/_230-270.sam} >> ${f/_filtered.sam/_230-270.header}; rm ${f/_filtered.sam/_230-270.sam}; mv ${f/_filtered.sam/_230-270.header} ${f/_filtered.sam/_230-270.sam}; samtools view -@ 23 -S -t /ix1/shainer/mouse/mm10/mm10.chrom.sizes -b -o ${f/_filtered.sam/_230-270.bam} ${f/_filtered.sam/_230-270.sam}; rm ${f/_filtered.sam/_230-270.sam}; done

for f in *_filtered.sam; do awk ' $9 <= 165 && $9 >= 135 || $9 >= -165 && $9 <= -135 ' $f > ${f/_filtered/_135-165}; cp bowtie2_mm10.header ${f/_filtered.sam/_135-165.header}; cat ${f/_filtered.sam/_135-165.sam} >> ${f/_filtered.sam/_135-165.header}; rm ${f/_filtered.sam/_135-165.sam}; mv ${f/_filtered.sam/_135-165.header} ${f/_filtered.sam/_135-165.sam}; samtools view -@ 23 -S -t /ix1/shainer/mouse/mm10/mm10.chrom.sizes -b -o ${f/_filtered.sam/_135-165.bam} ${f/_filtered.sam/_135-165.sam}; rm ${f/_filtered.sam/_135-165.sam}; done

for f in *_filtered.sam; do awk ' $9 <= 130 && $9 >= 100 || $9 >= -130 && $9 <= -100 ' $f > ${f/_filtered/_100-130}; cp bowtie2_mm10.header ${f/_filtered.sam/_100-130.header}; cat ${f/_filtered.sam/_100-130.sam} >> ${f/_filtered.sam/_100-130.header}; rm ${f/_filtered.sam/_100-130.sam}; mv ${f/_filtered.sam/_100-130.header} ${f/_filtered.sam/_100-130.sam}; samtools view -@ 23 -S -t /ix1/shainer/mouse/mm10/mm10.chrom.sizes -b -o ${f/_filtered.sam/_100-130.bam} ${f/_filtered.sam/_100-130.sam}; rm ${f/_filtered.sam/_100-130.sam}; done

for f in *_filtered.sam; do awk ' $9 <= 80 && $9 >= 1 || $9 >= -80 && $9 <= -1 ' $f > ${f/_filtered/_1-80}; cp bowtie2_mm10.header ${f/_filtered.sam/_1-80.header}; cat ${f/_filtered.sam/_1-80.sam} >> ${f/_filtered.sam/_1-80.header}; rm ${f/_filtered.sam/_1-80.sam}; mv ${f/_filtered.sam/_1-80.header} ${f/_filtered.sam/_1-80.sam}; samtools view -@ 23 -S -t /ix1/shainer/mouse/mm10/mm10.chrom.sizes -b -o ${f/_filtered.sam/_1-80.bam} ${f/_filtered.sam/_1-80.sam}; rm ${f/_filtered.sam/_1-80.sam}; done

#Sort and index files in preparation for danpos2 analysis
#See MNase_danpos.bash script for follow-up scripts 

module purge
module load gcc/8.2.0 samtools/1.14

for f in *270.sam; do samtools sort -@ 23 -O BAM -o ${f/.sam/.bam} $f; samtools index -@ 23 ${f/.sam/.bam}; done
for f in *165.sam; do samtools sort -@ 23 -O BAM -o ${f/.sam/.bam} $f; samtools index -@ 23 ${f/.sam/.bam}; done
for f in *130.sam; do samtools sort -@ 23 -O BAM -o ${f/.sam/.bam} $f; samtools index -@ 23 ${f/.sam/.bam}; done
for f in *80.sam; do samtools sort -@ 23 -O BAM -o ${f/.sam/.bam} $f; samtools index -@ 23 ${f/.sam/.bam}; done
