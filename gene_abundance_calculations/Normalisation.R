#Gene Normalisation Script
#Jason Bosch
#Last updated: 23/07/2020
#
#This script is based on the tutorial by Max Ortiz. It is intended to determine the relative abundance of genes of interest
#compared to reference single copy genes. It works using the result of hmmsearch.
#
#Script usage:
#
#Rscript Normalisation.R hmmsearch_results filtered_reads_file key_file significance_level outputname
#
#hmmsearch results is the file outputted by the --tblout flag.
#filtered_reads_file is the combined forward and reverse reads filtered with whichever criteria you are using
#key_file is a quoted csv table with no header and three columns: sequence name as used in the hmm model, the average amino acid length given by the hmm model and E/C which designates that the gene is either (E)xperimental or a single-copy (C)ontrol
#significance_level is a number in the format 1e-6. Reads with an E-value greater than that will be discarded.
#output name is whatever you want the result file to be called

##################
##################

#Read in the files of interest
input <- commandArgs(trailingOnly = T)
target_reads <- read.table(input[1],stringsAsFactors = F)
sample_read_file <- input[2]
key <- read.table(input[3],stringsAsFactors = F, quote = '"',sep = ",")
siglvl <- as.double(input[4])
outputname <- input[5]

#Calculate the counts for each gene after filtering for significance level
target_reads_filtered <- subset(target_reads, target_reads$V5<=siglvl)
for (n in 1:nrow(key)) {
  key[n,4] <- nrow(subset(target_reads_filtered,target_reads_filtered$V3 == key[n,1]))
}

#Calculate the total number of reads in sample
#For sample 39 = 104363534
sample_reads <- as.numeric(system(paste("grep -c '>' ",sample_read_file,sep = ""),intern = T))
scaling_factor <- sample_reads/1000000

#Generate RPKM
for (n in 1:nrow(key)) {
  key[n,5] <- (key[n,4]/scaling_factor)/((key[n,2]*3)/1000)
}

# #There are alternative equations for RPKM but they all give the same answer.
# #V1 - https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
# (key[n,4]/scaling_factor)/((key[n,2]*3)/1000)
# #V2 - www.metagenomics.wiki/pdf/definition/rpkm-calculation
# key[n,4]/(((key[n,2]*3)/1000)*sample_reads/1000000)
# #V3 - Max's tutorial
# (key[n,4]/scaling_factor)*((1000000)/(key[n,2]*3)/1000)

#Get normalisation value from the control genes
normalisation_factor <- exp(mean(log(key[key[,3]=="C",5])))

#Normalise gene counts
key[,6] <- key[,5]/normalisation_factor

#Save results
colnames(key) <- c("Gene","Length (aa)","E/C",paste("Mapped_Reads_@_",siglvl,sep = ""),"RPKM","Normalised_gene")
write.table(key,paste(outputname,".csv",sep = ""),quote = T, row.names = F, sep = ",")
