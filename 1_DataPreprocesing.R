###### DataPreprocessing.R ############################################
### OUTPUT: otu.no0 : filtered and 0's replaced otu table. 
###         g_otu.no0 : filtered and 0's replaced genus level otu table.
###         SAMPLEDATA : labels for samples 
###         TAX : taxa for otus
#######################################################################
# ASLI BOYRAZ, Dec 2021
#######################################################################

#Read data
library(readr)
tax_otu <- read.table('data/processed_taxatable.txt', sep="\t",header = TRUE) 
tax<-tax_otu[,1:7]
otu<-tax_otu[,-(1:7)]
labels <- read.table('data/task-healthy-cirrhosis.txt', sep="\t",header = FALSE) 


## Create Phyloseq Object

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
library(phyloseq)

OTU<-otu_table(otu, taxa_are_rows=TRUE) # taxas x samples
TAX<-tax_table(as.matrix(tax))
SAMPLEDATA<-sample_data(as.data.frame(labels))
taxa_names(OTU)<-taxa_names(TAX)
sample_names(OTU)<-sample_names(SAMPLEDATA)
cir<- phyloseq(OTU,TAX,SAMPLEDATA)


#Filtering:Taxa that were not seen with more than 20 counts in at least 30% of samples are filtered.
cir2 <-  filter_taxa(cir, function(x) sum(x > 20) > (0.3*length(x)), TRUE)

rm(otu,tax,labels)


# Get filtered data from phyloseq object.
OTU<-t(as.data.frame(otu_table(cir2))) # samples x taxas
SAMPLEDATA$Label<-ifelse(grepl("Cirrhosis", SAMPLEDATA$V2), 1, 0)
names(SAMPLEDATA)<-c("sampleID","Status","Label")
TAX<-as.data.frame(tax_table(cir))
TAX$otuid<-rownames(TAX)


#Replace 0's and close data.
library(zCompositions)
otu.no0<-cmultRepl(OTU) # samples x taxas
#rowSums(otu.no0) == 1

########################################
######   CREATE GENUS LEVEL OTU TABLE ##
########################################
install_github("microbiome/microbiome")
library(microbiome)
cir.genus <- aggregate_taxa(cir2, 'g')
g_OTU<-t(as.data.frame(otu_table(cir.genus))) #taxa x samples
g_otu.no0<-cmultRepl(g_OTU)
rm(g_OTU)

