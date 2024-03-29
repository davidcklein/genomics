# genomics
Master reference scripts, generalizable to many NGS-based analyses with minor modifications

# Parent directory setup
Each of the scripts housed in this repository assumes that the analysis will be run from a master directory ($parent), containing gzipped .fastq files in a subdirectory ($parent/raw). 

Over the course of the analysis, new subdirectories will be created within $parent, including those to house mapped alignment files ($parent/bam), normalized genomic coverage signal files ($parent/bw), and QC information ($parent/QC).

Because peak-calling and plotting can require much more specifically targeted analyses, I am not including these specific commands in the analysis scripts included in this repository. As it is, these scripts will run through the basic, repeatable analyses that are generalizable to most similar experiments. 

# Other information
The bash scripts are intended to be run in a SLURM environment, but are adaptable to other workload managers/cloud computing environments with adjustments to the header. 
R scripts are amenable to running in whatever R environment is preferred. It is assumed that all packages that are required for analysis are installed. Typically, these packages can be installed via CRAN or BioConductor repositories. 
