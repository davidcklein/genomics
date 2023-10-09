#H/T Ben Patty for this script
setwd("~/working_directory")

#Set up a table containing your counts as an R object
cts <- read.table("UN-NORMALIZED_counts_file.txt",header=TRUE,row.names=1)
head(cts)
coldata<-read.table("coldata.txt",header=TRUE)
head(coldata)
head(cts)

#Set up a TPM function to run the TPM normalization between the genes and their lengths

tpm <- function(cts, lengths) {
  rate <- cts / lengths
  rate / sum(rate) * 1e6
}

#Write in gene counts data as a df and apply TPM function on a gene-by-gene basis
genes_counts_Data <- data.frame(gene_name = row.names(cts), feature_length = as.integer(cts[,1]))
tpms <- apply(cts, 2, function(x) tpm(x, genes_counts_Data$feature_length))
head(tpms)
#Write a CSV file containing the gene identifiers and their TPM values
write.csv(as.data.frame(tpms), 
          file="Gene_counts_TPM.csv")

