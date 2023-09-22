#! /bin/bash
#SBATCH -N 1
#SBATCH -J RNAseq
#SBATCH -t 00-09:59
#SBATCH -c 4
#SBATCH --mem 100G
#SBATCH -o RNAseq_pipeline.out.txt
#SBATCH -p scavenger
##################################

module load gcc/8.2.0 star/2.7.0e samtools/1.14

hg38=/ix1/shainer/Dave/Master_Reference_Files/programs/STAR/hg38/
mm10=/ix1/shainer/Dave/Master_Reference_Files/programs/STAR/mm10/mm10_gencode


for f in *R1.fastq; do STAR --runThreadN 4 --genomeDir $hg38 --readFilesIn $f ${f/R1/R2} --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverReadLmax 0.02 --outFilterMultimapNmax 1 --outFileNamePrefix ${f/_R1.fastq/}; done

for f in *.bam; do mv $f ${f/Aligned.SortedByCoord/}; done

module purge
module load gcc/8.2.0 samtools/1.14

for f in *.bam; do samtools view -@ 3 -h -O BAM -q 7 -o ${f/.bam/_qf.bam} $f; samtools sort -@ 3 -o ${f/.bam/_sorted.bam} ${f/.bam/_qf.bam}; samtools index -@ 3 ${f/.bam/_sorted.bam}; done

module purge
module load subread/2.0.1

featureCounts -t transcript -s 2 -p -B -T 4 -a /ix1/shainer/Dave/Master_Reference_Files/hg38/gtf/hg38.gencode.v38.gtf -o Gencode_hg38_transcript_counts.txt *_sorted.bam
featureCounts -t gene -s 2 -p -B -T 4 -a /ix1/shainer/Dave/Master_Reference_Files/hg38/gtf/hg38.gencode.v38.gtf -g gene_name -o Gencode_hg38_gene_counts.txt *_sorted.bam

module purge
module load deeptools

for f in *_sorted.bam; do bamCoverage -b $f -o ${f/_sorted.bam/_F.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM --filterRNAstrand forward; done
for f in *_sorted.bam; do bamCoverage -b $f -o ${f/_sorted.bam/_R.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM --filterRNAstrand reverse; done
for f in *_sorted.bam; do bamCoverage -b $f -o ${f/_sorted.bam/_unstranded.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM; done
