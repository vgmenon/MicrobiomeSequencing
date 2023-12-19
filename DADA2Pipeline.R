#DADA2 Processing - Sample inference from amplicon data 
#Use to process the fastq files as provided by illumina sequencing

library(dada2); packageVersion("dada2")
library("knitr")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(DECIPHER); packageVersion("DECIPHER")
library(gridExtra); packageVersion("gridExtra")
library(phangorn); packageVersion("phangorn")
library(tidyverse); packageVersion("tidyverse") 
library(breakaway); packageVersion("breakaway")  
library(rms); packageVersion("rms")
library(dendextend); packageVersion("dendextend")   
library(HMP); packageVersion("HMP")       
library(picante); packageVersion("picante")  


#assign the working directory with the specific raw reads (Illumina 16S fastqz.gz files)
path <- "" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# Currently adapted to the Illumina reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Reads' quality profile
# Green: mean quality score
# Orange: quality quartiles
# Red: scaled proportion reads - not important for illumina reads
# Gray-scale: heat map of frequency of each quality score at each base position

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# forward reads are generally better quality - trim the last few nucleotides to control errors
# reverse reads are worse quality - common for illumina sequencing
# DADA2 incorporate the quality information - robust to lower quality sequencing.

# For V3-V4, trunclen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,270),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Error rates: parametric error model using every amplicon dataset - learnErrors method learns error model from data, 
# by alternatibe estimation of error rates and inference of sample composition until they converge on a consistent solution.

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# produces a plot for error rates for each possible nucleotide transition.
# each point is the observed error rate for each consensus quality score.
# black line - estimated error rates after convergence of the ML algorithm
# Red line: error rates expected under nominal definition of Q score
#Tldr: black line would ideally be a good fit to the points and error rates drop with quality
plotErrors(errF, nominalQ=TRUE)

# Sample inference: https://www.nature.com/articles/nmeth.3869#methods
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #Forward
dadaRs <- dada(filtRs, err=errR, multithread=TRUE) #Reverse
dadaFs[[1]] 

# Merge forward and reverse reads to obtain denoised full sequences
# align denoised forward reads with reverse compliment and constructing 'contig' sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[2]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

silva <- "" # assign the path for silva_nr99_v138.1_train_set.fa file for taxa assignment
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#making standard tables
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
ASV_seq <- "" # filename for the output of final ASV seqs
write(asv_fasta, ASV_seq)

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
ASV_table <- "" # output file name for the counts for each ASV
write.table(asv_tab, ASV_table, sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
asv_tax <- taxa
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

ASV_taxatable <- "" # output file for the ASV taxa table
write.table(asv_tax, ASV_taxatable, sep = "\t", quote=F, col.names=NA)









