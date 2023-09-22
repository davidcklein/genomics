#!/bin/bash
#
#SBATCH -N 1
#SBATCH --job-name=MNase
#SBATCH --output=MNase_Analysis_danpos_%A.txt
#SBATCH --time=05-23:59
#SBATCH --cpus-per-task=1
#SBATCH --mem=100g
#SBATCH -p scavenger
#############################

##Danpos2 is ANCIENT -- note that loading current modules (e.g. r/4.X.X or samtools/1.X) will break the software
#Danpos does not support multithreading in any way that I can tell, so it's not worth running this section with a lot of cores (hence the separate scripts for step 1 and step 2)
#Before running, change the sample names within the scripts to match the samples in your directory. I recommend keeping the same file name for all files, with the size class included at the end.

#set variable for danpos2 ($danpos2)
danpos2=/ihome/crc/install/danpos2/danpos-2.2.2/

module purge
module load gcc/8.2.0 r/3.5.1 samtools/0.1.19 danpos2

mkdir -p Danpos_1-80
mkdir -p Danpos_100-130
mkdir -p Danpos_135-165
mkdir -p Danpos_230-270

#Include all of your bam files, separated by commas, within the same command immediately following the dpos command
# -m 1 denotes a paired-end sequencing dataset
# -f 1 denotes that we will be using FDR values. If you do not want to use FDR values (will speed up the command but is not as good for calling nucleosome positions), set to 0
# -s 1 denotes that intermediate files will be saved. Set to 0 if you do not want these files saved. 
# --mifrsz and --mafrsz set min and max fragment sizes. Make sure that your min and max fragment sizes include the size class specified by your sam file -- including more fragment sizes is okay, since the awk command in the prior script will have size-selected already. 
# --extend will expand the sequencing reads to a length specified, similar to -e in deeptools

python $danpos2/danpos.py dpos SAMPLE_1_1-80.bam,SAMPLE_2_1-80.bam -m 1 -o Danpos_1-80 -f 1 -s 1 -m 1 --mifrsz 0 --mafrsz 100 --extend 80
python $danpos2/danpos.py dpos SAMPLE_1_100-130.bam,SAMPLE_2_100-130.bam -m 1 -o Danpos_100-130 -f 1 -s 1 -m 1 --mifrsz 100 --mafrsz 200 --extend 80
python $danpos2/danpos.py dpos SAMPLE_1_135-165.bam,SAMPLE_2_135-165.bam-m 1 -o Danpos_135-165 -f 1 -s 1 -m 1 --mifrsz 100 --mafrsz 200 --extend 80
python $danpos2/danpos.py dpos SAMPLE_1_230-270.bam,SAMPLE_2_230-270.bam -m 1 -o Danpos_230-270 -f 1 -s 1 -m 1 --mifrsz 200 --mafrsz 300 --extend 80

mkdir -p Danpos_1-80/wiq/
mkdir -p Danpos_100-130/wiq/
mkdir -p Danpos_135-165/wiq/
mkdir -p Danpos_135-165/wiq/

#Wiq is a quantile-normalization command intended to normalize five quintiles of data across multiple datasets within a single set of experiments. This can reduce variability in datasets and ensure that highly-occupied regions are not altering normalization of lowly occupied regions between samples. 

#buffer_size specifies the memory limit. Do not exceed ~80% of the memory requested in the slurm header. 

python $danpos2/danpos.py wiq /ix1/shainer/mouse/mm10/mm10.chrom.sizes ./Danpos_1-80/pooled/ --out_dir ./Danpos_1-80/wiq/ --buffer_size 80
python $danpos2/danpos.py wiq /ix1/shainer/mouse/mm10/mm10.chrom.sizes ./Danpos_100-130/pooled/ --out_dir ./Danpos_100-130/wiq/ --buffer_size 80
python $danpos2/danpos.py wiq /ix1/shainer/mouse/mm10/mm10.chrom.sizes ./Danpos_135-165/pooled/ --out_dir ./Danpos_135-165/wiq/ --buffer_size 80
python $danpos2/danpos.py wiq /ix1/shainer/mouse/mm10/mm10.chrom.sizes ./Danpos_230-270/pooled/ --out_dir ./Danpos_230-270/wiq/ --buffer_size 80

#Call nucleosome positions on quantile-normalized data:
python $danpos2/danpos.py dpos ./Danpos_1-80/wiq/ -o ./Danpos_1-80/qnor_result/ -f 1 --mifrsz 0 --mafrsz 100 --extend 80
python $danpos2/danpos.py dpos ./Danpos_100-130/wiq/ -o ./Danpos_100-130/qnor_result/ -f 1 --mifrsz 100 --mafrsz 200 --extend 80
python $danpos2/danpos.py dpos ./Danpos_135-165/wiq/ -o ./Danpos_135-165/qnor_result/ -f 1 --mifrsz 100 --mafrsz 200 --extend 80
python $danpos2/danpos.py dpos ./Danpos_230-270/wiq/ -o ./Danpos_230-270/qnor_result/ -f 1 --mifrsz 200 --mafrsz 300 --extend 80
