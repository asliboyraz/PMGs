# PMGs

*Principal Microbial Groups : a new way to explore microbiome data*

***
PMG = Principal Microbial Groups.
This R package provides functions for the analysis of compositional data (e.g., data representing proportions of different variables/parts). Specifically, this package allows analysis of microbiome data where the OTUs can be grouped through a SBP (sequential binary partitioning) and provides a balance of OTU groups to be utilized for for the search of biomarkers in human microbiota. 

## Authors ##
Aslı Boyraz, Middle East Technical University

## Quick Start ##

`PMGs` procedure is illustrated on a cirrhosis study data. Data is obtained from the https://github.com/knights-lab/MLRepo.

Load Cirrhosis Data
``` r
library(readr)
tax_otu <- read.table('data/processed_taxatable.txt', sep="\t",header = TRUE) 
tax<-tax_otu[,1:7]
otu<-tax_otu[,-(1:7)]
labels <- read.table('data/task-healthy-cirrhosis.txt', sep="\t",header = FALSE) 
```

A Phyloseq object created. Data is filtered and 0's are replaced. TAX table and Sampledata are organized.

```{r}
library(phyloseq)
OTU<-otu_table(otu, taxa_are_rows=TRUE) # taxas x samples
TAX<-tax_table(as.matrix(tax))
SAMPLEDATA<-sample_data(as.data.frame(labels))
taxa_names(OTU)<-taxa_names(TAX)
sample_names(OTU)<-sample_names(SAMPLEDATA)
cir<- phyloseq(OTU,TAX,SAMPLEDATA)
#Taxa that were not seen with more than 20 counts in at least 30% of samples are filtered.
cir2 <-  filter_taxa(cir, function(x) sum(x > 20) > (0.3*length(x)), TRUE)
OTU<-t(as.data.frame(otu_table(cir2)))
library(zCompositions)
otu.no0<-cmultRepl(OTU)
TAX<-as.data.frame(tax_table(cir))
TAX$otuid<-rownames(TAX)
SAMPLEDATA$Label<-ifelse(grepl("Cirrhosis", SAMPLEDATA$V2), 1, 0)
names(SAMPLEDATA)<-c("sampleID","Status","Label")
```

Before `PMGs` construction, the analyst should decide the minimum number of groups to construct.

```{r}
V<- mPBclustvar(X) #
coord<-milr(X,V)
plot(diag(var(coord)),xlab="number of coordinates",ylab="Explained Variance")
```
![](README-plot-1.png)











## Bugs/Feature requests ##
I appreciate bug reports and feature requests. Please post to the github issue tracker [here](https://github.com/asliboyraz/pmgs/issues). 


