
library("phyloseq"); packageVersion("phyloseq") #1.38.0
library("ggplot2"); packageVersion("ggplot2") #3.3.5
library(readr)
#install.packages("lubridate")
library(lubridate)
library("plyr")
#install.packages("dplyr")
library(dplyr)
library(ggpubr)
library(tibble)
library(vegan)
library("vegan")
library("tidyr")
library("reshape")
library("gplots")
#install.packages("openxlsx")
library(openxlsx)
#install.packages("viridis")
library(viridis)
#install.packages("readxl")
library(readxl)
#install.packages("writexl")
library(writexl)
library(dada2); packageVersion("dada2")
library("knitr")
library(Biostrings); packageVersion("Biostrings")
library(DECIPHER); packageVersion("DECIPHER")
library(gridExtra); packageVersion("gridExtra")
library(phangorn); packageVersion("phangorn")
library(tidyverse); packageVersion("tidyverse") 
library(DESeq2); packageVersion("DESeq2")         
library(breakaway); packageVersion("breakaway")  
library(rms); packageVersion("rms")
library(dendextend); packageVersion("dendextend")   
library(HMP); packageVersion("HMP")      
library(picante); packageVersion("picante") 
library(RColorBrewer)
install.packages("remotes")
library(remotes)
devtools::install_github("vmikk/metagMisc")

# Reading the ASV table - produced by running the fastq illumina raw reads through the DADA2 pipeline
ASV_table <- "" # type in the .tsv file name obtained from DADA2 sequencing
count_tab_VGM <- read.table(ASV_table, header=T, row.names=1, check.names=F, sep="\t")

ASV_counts_VGM <- as.matrix(count_tab_VGM) #convert data table into matrix
ASV_counts_VGM.df <- as.data.frame(ASV_counts_VGM) #convert matrix into data frame

# matrix has a single class, whereas dataframe can have multiple classes of data


# Assign a taxa to each ASV - ASV taxa file obtained from DADA2 pipeline
ASV_taxa <- "" #input the file name for the ASV taxa file obtained from DADA2 pipeline
asv_tax_VGM <- read.delim(ASV_taxa, header=FALSE, fill=TRUE, row.names = 1)

asv_tax_VGM.df <- asv_tax_VGM

#Assign column names as the titles in row 1
colnames(asv_tax_VGM.df) <- asv_tax_VGM[1, ]  
#Drop row 1
asv_tax_VGM.df <- asv_tax_VGM.df[-1,]

# results in a ASV taxa file with column names as the taxonomical rank.

asv_taxm_VGM <- as.matrix(asv_tax_VGM.df)

#make ASV_counts numeric
ASV_countsn_VGM <- lapply(ASV_counts_VGM.df, as.numeric)

# obtain the metadata for the sample treatments 
Meta <- "" # input the file name for the csv file containing the metadata for the samples sequenced
VGM_metadata <- read_csv(Meta)
# convert row names to sample names
VGM_metadata <- VGM_metadata %>% remove_rownames %>% column_to_rownames(var="Sample")

# create a phyloseq object from the ASV counts, metadata and ASVC taxa map.
VGM_phyloseq <- phyloseq(otu_table(ASV_counts_VGM.df, taxa_are_rows=TRUE), 
                         sample_names(VGM_metadata), 
                         tax_table(asv_taxm_VGM))

#check the phyloseq object to make sure your samples align

#convert everything to % relative abundance
VGM_phyloseq_rel = transform_sample_counts(VGM_phyloseq, function(x) (x / sum(x))*100 )

############################### Level ###################################
# Modify the taxa level as desired                                        
taxa = "Class" 
taxan = 3
 
taxa_counts_tab_VGM <- otu_table(tax_glom(VGM_phyloseq, taxa))
head(taxa_counts_tab_VGM)
# making a vector of Phylum names to set as row names
phyla_tax_vec_VGM <- as.vector(tax_table(tax_glom(VGM_phyloseq, taxrank=taxa))[,taxan])
rownames(taxa_counts_tab_VGM) <- as.vector(phyla_tax_vec_VGM)

dim(taxa_counts_tab_VGM)

taxa_counts_tab_VGM <- transform_sample_counts(taxa_counts_tab_VGM, function(x) (x / sum(x))*100 )
  
#make a different file for manipulating
major_taxa_for_plot_VGM <- as.data.frame(taxa_counts_tab_VGM)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot_VGM$MTaxa <- row.names(major_taxa_for_plot_VGM)

#gather function converts multiple columns into a single column key-pair
major_taxa_for_plot_VGM.g <- gather(major_taxa_for_plot_VGM, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot_VGM.g)

# now we want a table with metadata for each sample
sample_info_for_merge_VGM<-data.frame("ID"=row.names(VGM_metadata),
                                      "fraction"=VGM_metadata$SizeFrac,
                                      stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
# add metadata to the major taxa table
major_taxa_for_plot_VGM.g2 <- merge(major_taxa_for_plot_VGM.g, sample_info_for_merge_VGM)
dim(major_taxa_for_plot_VGM.g2)

outname <- "" # input the file name for the output of major taxa
write.xlsx(major_taxa_for_plot_VGM, paste0(outname,taxa,'.xlsx'))

phylumGlommed = tax_glom(VGM_phyloseq_rel, taxa)

p <- plot_bar(phylumGlommed, fill =taxa)+geom_bar(stat="identity")#+xlab(">3um Size Fraction")+ylab("Relative Abundance [%]")
print(p)
############################### Filtering ######################################

###############################Genera Level ####################################

new_gplot_VGM <- major_taxa_for_plot_VGM.g2[major_taxa_for_plot_VGM.g2$MTaxa =="Proteobacteria" |
                                              major_taxa_for_plot_VGM.g2$MTaxa =="Bacteroidota" |
                                              major_taxa_for_plot_VGM.g2$MTaxa =="Actinobacteriota" |
                                              major_taxa_for_plot_VGM.g2$MTaxa =="Cyanobacteria" |
                                              major_taxa_for_plot_VGM.g2$MTaxa =="Firmicutes" |
                                              major_taxa_for_plot_VGM.g2$MTaxa =="Marinimicrobia (SAR406 clade)",]



write.xlsx(new_gplot_VGM, 'Relativeabundance_MajorPhylum.xlsx')

############################### Family Level ####################################

new_gplot_family_VGM <- major_taxa_for_plot_VGM.g2[major_taxa_for_plot_VGM.g2$MTaxa =="Rhodobacteraceae" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa =="Colwelliaceae" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa =="Clade I" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa =="Actinomarinaceae" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "SAR86 clade" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Flavobacteriaceae" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Clade II"|
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Oxalobacteraceae"|
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Alteromonadaceae"|
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Porticoccaceae"|
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Enterobacteriaceae",]

write.xlsx(new_gplot_family_VGM, 'Relativeabundance_MajorFamily.xlsx')

###############################Genera Level ####################################

new_gplot_genus_VGM <- major_taxa_for_plot_VGM.g2[major_taxa_for_plot_VGM.g2$MTaxa =="Amylibacter" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa =="Colwellia" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa =="Cognatiyoonia" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa =="Alloprevotella" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Aurantivirga" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Alteromonas" |
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Clade Ia"|
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Aurantivirga"|
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Ulvibacter"|
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Sulfitobacter"|
                                                     major_taxa_for_plot_VGM.g2$MTaxa == "Pseudofulvibacter",]

write.xlsx(new_gplot_genus_VGM, 'Relativeabundance_MajorGenus.xlsx')

