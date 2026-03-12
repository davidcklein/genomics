#! /bin/bash
#SBATCH -N 1
#SBATCH -J RNAseq
#SBATCH -t 00-09:59
#SBATCH -c 4
#SBATCH --mem 48G
#SBATCH -o RNAseq_pipeline_%A.out.txt
##################################
module load gcc star/2.7.0e samtools

hg38=/ix1/shainer/Dave/Master_Reference_Files/programs/STAR/hg38/
mm10=/ix1/shainer/Dave/Master_Reference_Files/programs/STAR/mm10/mm10_gencode
T2T=/ix1/shainer/Dave/Master_Reference_Files/programs/STAR/T2T

for f in *R1.fastq.gz; do if [[ ! -e "${f/_R1.fastq.gz/_Aligned.sortedByCoord.out.bam}" ]]; then STAR --runThreadN 4 --genomeDir $T2T --readFilesIn $f ${f/R1/R2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverReadLmax 0.1 --outFilterMultimapNmax 1 --outFileNamePrefix ${f/_R1.fastq.gz/_}; fi; done
#Align read fragments using a splice-aware aligner (I prefer STAR, but Hisat2 is another good option. I would avoid TopHat altogether at this point).

for f in *Aligned.sortedByCoord.out.bam; do mv $f ${f/Aligned.sortedByCoord.out/}; done

for f in *.bam; do if [[ ! -e "${f/.bam/_qf.bam}" ]]; then samtools view -@ 3 -h -O BAM -q 7 -o ${f/.bam/_qf.bam} $f; samtools index -@ 3 ${f/.bam/_qf.bam}; fi; done
#This command will filter out reads with a MAPQ alignment quality score of < 7 and index the resulting bam files

rm *.bam *.bai
for f in *_qf.bam; do mv $f ${f/_qf/}; done
for f in *.bam; do samtools index -@ 3 $f; done

module purge
module load subread
if [[ ! -e "NCBI_Refseq_T2T_gene_counts.txt" ]]; then featureCounts -t gene -s 2 -p --countReadPairs -B -T 4 -a /ix1/shainer/human/chm13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf -g gene_id -o NCBI_Refseq_T2T_gene_counts.txt *.bam; fi
#This command will use subread to assign reads to genes for T2T.

module purge
module load deeptools
for f in *.bam; do
    if [[ ! -e "${f/.bam/_F.bw}" ]]; then bamCoverage -b $f -o ${f/.bam/_F.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM --filterRNAstrand forward; fi
    if [[ ! -e "${f/.bam/_R.bw}" ]]; then bamCoverage -b $f -o ${f/.bam/_R.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM --filterRNAstrand reverse; fi
    if [[ ! -e "${f/.bam/_unstranded.bw}" ]]; then bamCoverage -b $f -o ${f/.bam/_unstranded.bw} -of bigwig -bs 5 --smoothLength 20 -p 4 --normalizeUsing BPM; fi
done
#These three deeptools commands will generate TPM (BPM) normalized strand-specific and unstranded genome coverage files. Each file is set to be binned in 5-bp increments and smoothed by averaging 20-bp regions. For higher resolution, decrease -bs to 1; for faster processing, increase the -bs parameter. Maintain a 4:1 ratio between smoothLength and bin size.
