#! /bin/bash
#SBATCH -N 1
#SBATCH -J RNAseq
#SBATCH -t 00-09:59
#SBATCH -c 4
#SBATCH --mem 100G
#SBATCH -o RNAseq_pipeline_%A.out.txt
##################################

module load gcc star/2.7.0e samtools

hg38=/ix1/shainer/Dave/Master_Reference_Files/programs/STAR/hg38/
mm10=/ix1/shainer/Dave/Master_Reference_Files/programs/STAR/mm10/mm10_gencode
T2T=/ix1/shainer/Dave/Master_Reference_Files/programs/STAR/T2T

for f in *R1.fastq; do STAR --runThreadN 4 --genomeDir $T2T --readFilesIn $f ${f/R1/R2} --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverReadLmax 0.02 --outFilterMultimapNmax 1 --outFileNamePrefix ${f/_R1.fastq/_}; done
#Align read fragments using a splice-aware aligner (I prefer STAR, but Hisat2 is another good option. I would avoid TopHat altogether at this point). 

for f in *.bam; do mv $f ${f/Aligned.sortedByCoord.out/}; done

module purge
module load gcc samtools

for f in *.bam; do samtools view -@ 3 -h -O BAM -q 7 -o ${f/.bam/_qf.bam} $f; samtools sort -@ 3 -o ${f/.bam/_sorted.bam} ${f/.bam/_qf.bam}; samtools index -@ 3 ${f/.bam/_sorted.bam}; done
#This command will filter out reads with a MAPQ alignment quality score of < 7, as well as coordinate-sort and index the resulting bam files

module purge
module load subread

featureCounts -t gene -s 2 -p --countReadPairs -B -T 6 -a /ix1/shainer/human/chm13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf -g gene_id -o NCBI_Refseq_T2T_gene_counts.txt *_sorted.bam

#This command will use subread to assign reads to genes for T2T.

module purge
module load deeptools

for f in *_sorted.bam; do bamCoverage -b $f -o ${f/_sorted.bam/_F.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM --filterRNAstrand forward; done
for f in *_sorted.bam; do bamCoverage -b $f -o ${f/_sorted.bam/_R.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM --filterRNAstrand reverse; done
for f in *_sorted.bam; do bamCoverage -b $f -o ${f/_sorted.bam/_unstranded.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM; done

#These three deeptools commands will generate TPM (BPM) normalized strand-specific and unstranded genome coverage files. Each file is set to be binned in 5-bp increments and smoothed by averaging 20-bp regions. For higher resolution, decrease -bs to 1; for faster processing, increase the -bs parameter. Maintain a 4:1 ratio between smoothLength and bin size. 
