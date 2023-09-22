# genomics
Master reference scripts, generalizable to many NGS-based analyses with minor modifications


# Parent directory setup
Each of the .bash scripts housed in this repository assumes that the analysis will be run from a master directory ($parent), containing gzipped .fastq files in a subdirectory ($parent/raw). 

In the course of the analysis, new subdirectories will be created within $parent, including those to house mapped alignment files ($parent/bam), normalized signal files ($parent/bw), peak files ($parent/peaks), QC information ($parent/QC), and heatmaps anchored over: RefSeq Select mRNA TSSs and gene-distal DNaseI hypersensitive sites (inferred to be likely enhancer regions). 


# Other information
The scripts are intended to be run in a SLURM environment, but are adaptable to other workload managers/cloud computing environments with adjustments to the header. 
