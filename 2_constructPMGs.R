############################################################################
###### constructPMGs.R #####################################################
### OUTPUT  PMGs table : Each column represents a PMG. 
###         OGT table :   (OtuID|Group|Taxa) information	
###         CODA Dendrogram on PMGs is available
#############################################################################
# NOTE 1: Call DataProcessing.R to obtain otu table without 0s. (otu.no0)
# NOTE 2: Before running findOptimalNumOfGroups method
#         - decide minNumberOfGroups looking at scree plot.
#         - maximumNumberOfGroups is handled in findOptimalNumOfGroups method 
#           and finds optimal number for groups.
##############################################################################
# ASLI BOYRAZ, Dec 2021
##############################################################################

source("lib/PMG_lib.R")

X<-otu.no0
V<- mPBclustvar(X) #
coord<-milr(X,V)
plot(diag(var(coord)),xlab="number of coordinates",ylab="Explained Variance") #scree plot

X<-otu.no0
#minimumNumOfGroups<- 70
#numberofgroup<-findOptimalNumOfGroups(otu.no0,minimumNumOfGroups)
numberofgroup<-27 # this is optimal for cirrhosis dataset.

##########################################################################
#### createPMGs : outputs PMGs table and OGT (OtuID|Group|Taxa) dataframe
##########################################################################
PMGs<-createPMGs(otu.no0,numberofgroup)  
rm(X,V,coord)

#########################################################################
#### Draw CODA Dendrogram on PMGs
#########################################################################
PMGs_labels<-cbind(PMGs,SAMPLEDATA$Status)
names(PMGs_labels)[ncol(PMGs_labels)]<-"label"
PMGs_labels$label<-as.factor(PMGs_labels$label)
W<- mPBclustvar(PMGs)
library(compositions)
CoDaDendrogram(X=acomp(PMGs),V=W,type="lines",range=c(-10,10))
CoDaDendrogram(X=acomp(PMGs[PMGs_labels$label=="Cirrhosis",]), col="red",add=TRUE,V=E,type="lines",range=c(-7,7))
CoDaDendrogram(X=acomp(PMGs[PMGs_labels$label=="Healthy",]), col="green",add=TRUE,V=E,type="lines",range=c(-7,7))
rm(PMGs_labels,W)
