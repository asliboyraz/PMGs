#########################################################################################################################
###############   CONSTRUCTING PRINCIPAL MICROBIAL GROUPS FUNCTIONS  ####################################################
#########################################################################################################################
# Functions
# 1- findOptimalNumOfGroups(X) : returns optimal number of groups to construct based on otutable(X). x must be non-zero.
# 2- findMaxNumOfGroups(X) : returns max number of groups for otutable(X). 
#                            It is decided based on the mean of coordinate variances.
#                            X must be non-zero.
# 3- createPMGs(X,numberofgroup) : returns PMG matrix.x must be non-zero.
# 4- groupOTUs(X,numberofgroup) : Groups otus and returns a table of otuid, taxa and group info.
#                                 x must be non-zero.
# 5- processColumns(x,labels) : Process columns on Contrast Matrix of Principle Balances(SignMatrix) for grouping.
#                               Returns signMatrix with group column added.
# 6- CheckPrevColsPOSorNEG(M,colno) : Used in processedColumns to process each column.
#
#########################################################################################################################

###################################################################
###### function findOptimalNumOfGroups ############################
### Finds Optimal Number of Principal Microbial Groups to construct 
#-----
# X: otu table without 0's.
# min: minimum number of groups decided by user.
#-----
# used in : -
# calls: findMaxNumOfGroups, createPMGs
###################################################################
# ASLI BOYRAZ, Dec 2021
###################################################################
findOptimalNumOfGroups <- function(X,min){
  library(dplyr) 
  library(stringr) 
  library(caTools) 
  library(caret)
  
  maxnumberofgroups <- findMaxNumOfGroups(X)
  minnumberofgroups <- min
  
  fit.control <- trainControl(method = "repeatedcv", number = 10, repeats = 10,
                              summaryFunction = twoClassSummary, classProbs = TRUE, allowParallel = F)
  
  roc.valuesMin<-NULL
  
  for(numberofgroup in c(minnumberofgroups:maxnumberofgroups)){
    
    matrixPMGs<-createPMGs(X,numberofgroup)
    
    V<- mPBclustvar(matrixPMGs)
    coord<-milr(matrixPMGs,V)
    
    mydataPGMs.ILR<-cbind(as.data.frame(coord),SAMPLEDATA$Status)
    names(mydataPGMs.ILR)[ncol(mydataPGMs.ILR)]<-"label"
    mydataPGMs.ILR$label<-as.factor(mydataPGMs.ILR$label)
    set.seed(123)
    fit.PMGbalance <- train(label ~ ., data = mydataPGMs.ILR, method = "glm", 
                            family = "binomial", trControl = fit.control)
    roc<-fit.PMGbalance$results[,"ROC"]
    roc.valuesMin <- c(roc.valuesMin,roc)
    print(numberofgroup)
  }
  
  roc.values<<-roc.valuesMin
  plot(roc.valuesMin,xlab=(c("Number of Groups")),ylab="AUC value", xaxt="n")
  axis(1, at=1:((maxnumberofgroups-minnumberofgroups)+1), labels=c(minnumberofgroups:maxnumberofgroups))
 
  optimalNumberofGroup<-as.numeric(which(roc.values %in% max(roc.values))) + (minnumberofgroups-1)
  
  print(paste0("Optimal Number of group is:",optimalNumberofGroup))
  return(optimalNumberofGroup)

}

#######################################################################
###### function findMaxNumOfGroups ####################################
### Finds maximum number of Principal Microbial Groups 
### to construct based on mean of coordinate variances 
#-----
# X: otu table without 0's.
#-----
# used in : -
# calls: mPBclustvar
#######################################################################
# ASLI BOYRAZ, Dec 2021
#######################################################################
findMaxNumOfGroups <- function(X){
  V<- mPBclustvar(X)
  coord<-milr(X,V)
  #mean(diag(var(coord)))
  ## number of coordinates that have higher value than mean value.
  ## Dont take account the points lower than mean value of var(coords).
  diag(var(coord))>mean(diag(var(coord)))
  maxnumberofgroups <- as.data.frame((table(diag(var(coord))>mean(diag(var(coord))))["TRUE"]))[1,1] + 1
  return(maxnumberofgroups)
}

###################################################################
###### function createPMGs ########################################
### Creates Principal Microbial Groups
#-----
# X: otu table without 0's.
# numberofgroup: number of Principal Microbial Groups to create
#-----
# used in : -
# calls: groupOTUs
###################################################################
# ASLI BOYRAZ, Dec 2021
###################################################################
createPMGs <- function(X,numberofgroup) {
  library(compositions)
  OGT<<-groupOTUs(X,numberofgroup)
  PMGs<-as.data.frame((matrix(ncol=numberofgroup,nrow=nrow(X))))  # 
  groupnames<-(paste0("G",1:numberofgroup))
  for(i in groupnames)
  {
    if(ncol(as.data.frame(X[,OGT$otuid[OGT['group'] == i]]))==1)
    {PMGs[[i]]<-X[,OGT$otuid[OGT['group'] == i]] }
    else{
      PMGs[[i]]<-geometricmeanRow(X[,OGT$otuid[OGT['group'] == i]]) }
  }
  PMGs[,c(1:numberofgroup)]<-NULL
  return(PMGs)
}

###########################################################################
###### function groupOTUs #########################################
### Creates groups out of otus 
#-----
# X: otu table without 0's.
# numberofgroup: number of Principal Microbial Groups to create
#-----
# used in : createPMGs
# calls: processColumns
###################################################################
# ASLI BOYRAZ, July 2018
###################################################################
groupOTUs <- function(X,numberofgroup) {
  groupnames<<-rev(paste0("G",1:numberofgroup))
  V<- mPBclustvar(X)
  coord<-milr(X,V)
  tempV<-V[,1:numberofgroup-1]
  s<-as.data.frame(sign(tempV))
  s$group<-NA
  s<-processColumns(s,numberofgroup)
  
  OTUsGroups<- cbind(rownames(s),s$group)
  colnames(OTUsGroups)<-c("otuid","group")
  OGT<-merge(OTUsGroups,TAX, by="otuid")
  return(OGT)
}

#####################################################################################
###### function processColumns ######################################################
### Process columns on Contrast Matrix of Principle Balances(SignMatrix) for grouping.
### Returns signMatrix with group column added.
#-----
# s: Contrast Matrix of Principle Balances (SignMatrix)
# numberofgroup: how many groups to construct out of otus.
#-----
# used in : groupOTUs
# calls: CheckPrevColsPOSorNEG0
####################################################################################
# ASLI BOYRAZ, July 2018
#####################################################################################

processColumns<-function(s,numberofgroup){
  for(i in rev(c(1:(numberofgroup-1)))){
    #print("---PROCESSING---")
    #print(i)
    
    if(i==(numberofgroup-1)){
      
      s[which(s[,i]<0),]$group <- groupnames[1]
      #print(groupnames[1])
      groupnames<<-groupnames[-1]
      # print("---group created---")
      s[which(s[,i]>0),]$group <- groupnames[1]
      #print(groupnames[1])
      groupnames<<-groupnames[-1]
      #print("---group created---")
    }
    
    else{
      flag<-CheckPrevColsPOSorNEG0(s,i)
      if(length(unique(flag[,1]))==1 && unique(as.list(flag[,1]))==TRUE){
        s[which(s[,i]>0),]$group <- groupnames[1]
        #print(groupnames[1])
        groupnames<<-groupnames[-1] 
      }
      if(length(unique(flag[,2]))==1 && unique(as.list(flag[,2]))==TRUE){
        
        s[which(s[,i]<0),]$group <- groupnames[1]
        #print(groupnames[1])
        groupnames<<-groupnames[-1] 
      }
    }
  }
  
  return(s)
}
#########################################################################
###### function CheckPrevColsPOSorNEG0 ##################################
### Check if previous columns for the rows with positive values and 
### negative values are all 0's.  
### Returns a flag with two columns: 
### -> First column is flags for rows with positive values. 
### -> Second column is flags for rows with negative values.
### Returns a TRUE column if previous columns only contain 0's.  
#-----
# M: Contrast Matrix of Principle Balances
# colno: the processed column no
#-----
# used in : processColumns
#########################################################################
# ASLI BOYRAZ, July 2018
#########################################################################

CheckPrevColsPOSorNEG0<-function(M,colno){
  pos_prevflag<-NULL
  neg_prevflag<-NULL
  
  POS<-which(M[,colno]>0)
  NEG<-which(M[,colno]<0)
  
  collist<-c((colno+1):(ncol(M)-1))
  
  for(i in collist) {
    #print(paste0("Checking for Col.",i))
    if (length(unique(M[POS,i])) < 2 )  { # POS rows on colno is All 0's
      pos_prevflag<- c(pos_prevflag,TRUE) #If all prev cols 0.
    }else{
      pos_prevflag<- c(pos_prevflag,FALSE) ##need to change column
    }
    
  }
  for(i in collist) {
    
    if (length(unique(M[NEG,i])) < 2 )  { # NEG rows on colno is All 0's
      neg_prevflag<- c(neg_prevflag,TRUE) #If all prev cols 0.
    }else{
      neg_prevflag<- c(neg_prevflag,FALSE) ##need to change column
    }
  }
  prevflag<-cbind(pos_prevflag,neg_prevflag)
  return(prevflag)
}
#####################################################################################
#####################################################################################
#####################################################################################
##################### FROM MINI COMPOSITION by Juanjo Egozcue #######################
#####################################################################################

###############################################
# function mPBclustvar(x)
# given a CoDa data set x[,1:D] a cluster analysis
# of the variables (columns) is carried out using
# Ward's method in stats{hclust}
# no zeros or negatives allowed
# The contrast matrix is computed and returned
#################################################
# JJ Egozcue, August 2017
#################################################
mPBclustvar<-function(x){
  # distances from variation
  hdist = sqrt(mvariation(x))
  hdist=as.dist(hdist)
  # cluster
  hmerge<-hclust(d=hdist,method="ward.D2")
  #plot(hmerge)
  # signs from $merge; contrast matrix
  Vsigns=mmerge2sign(hmerge$merge)
  V = mbuildcontrast(Vsigns)
  rownames(V)=colnames(x)
  return(V)
}
#################################################

###################################################
# function mmerge2sign transform the tree code $merge
# obtained in a hierachical cluster into a sign
# code of a SBP. The result can be used to construct
# a contrast matrix V using mbuildConstrast
# It is based on the function gsi.merge2signary in
# package compositions.
# Merge a (D,2) matrix of type $merge
# V on the output a SBP sign code (D,D-1) matrix
# containing the SBP code.
##############################################
# JJ Egozcue, August 2017
##############################################
mmerge2sign<-function(Merge){
  V = matrix(0, ncol = nrow(Merge) + 1, nrow = nrow(Merge))
  for (i in 1:nrow(Merge)) {
    for (j in 1:2) {
      weight = (-1)^j
      k = Merge[i, j]
      if (k < 0) {
        V[i, abs(k)] = weight
      }
      if (k > 0) {
        take = as.logical(V[k, ])
        V[i, take] = rep(weight, sum(take))
      }
    }
  }
  revV = V[nrow(V):1,]
  return(t(revV))
}
################################################

###### function mvariation #######
## computes the variation matrix
## of a compositional data matrix x
## There are three options of normalization
## optnor="none"  variation matrix as is
## optnor="minassoc" variation over minimum
##        associated matrix
## optnor=" ", rational 0,1 transform
## These normalizations are
## "minassoc" : t_ij = (d-1) v_ij / 2 totvar
## "linear01" : t_ij = totvar / [totvar + ((d-1)/2) v_ij]
## The value of v_ij which distributes uniformly total
## variance is v_u = 2 totvar/(d-1)
## characteristic values
##  "no normalization"   "minassoc"  "linear01"
##     0                    0             1
##     v_u                  1             0.5
##     Inf                  Inf           0
#########################################################
# JJ Egozcue, May 2017
###########################
mvariation <- function(x,optnor="none"){ 
  if(is.null(dim(x)) ){x=t(as.matrix(x))}
  d = ncol(x)
  xclr = mclr(x)
  co = var(xclr)
  va = diag(co)
  co1 = matrix(rep(va,each=d),ncol = d)
  co2 = matrix(rep(va,times=d),ncol = d)
  varia = -2 * co + co1 + co2
  if(optnor=="none"){  return(varia) }
  totvar = sum(varia)/(2*d)
  varian= varia*((d-1)/(2*totvar))
  if(optnor=="minassoc"){
    return(varian)
  }
  if(optnor=="linear01"){
    varian1 = totvar/(varia*((d-1)/2)+totvar)
    return(varian1)
  }
}


##### function mclr ##############
## computes the clr of x by rows
###########################
# JJ Egozcue, May 2017
###########################
mclr <- function(x){
  if(is.null(dim(x)) ){x=t(as.matrix(x))}
  logx = log(x)
  rs = outer(rowSums(logx),rep((1/ncol(x)),length=(ncol(x))))
  xclr=logx-rs
  return(xclr)
}

###### function milr ##############
## computes the ilr transformation
## of x  (n,D)
## using the contrast matrix V
###########################
# JJ Egozcue, May 2016
###########################
milr <- function(x, V){
  if(is.null(dim(x)) ){x=t(as.matrix(x))}
  xilr = log(as.matrix(x)) %*% V
  return(xilr)
}

####################################
##### function mbuildcontrast ###
##### Builds the contrast matrix ###
##### from a SBP code (signs) ######
##### for an ilr-transform    ######
# contrasts are defined by columns
# intput W is a matrix with D rows
#     and a undefined number of columns
#     containing 1,-1,0
# W also can be a vector e.g. c(1,1,-1,0)
# The matrix W does not need to correspond
# to a sequential binary partition,
# or to have strictly D-1 columns. Each
# row define a balance and it is computed
# independently of other columns.
####################################
## J. J. Egozcue, May 2014, May 2017
####################################
mbuildcontrast <-function(W=c(1,-1)){
  if(length(W)<2){
    print("improper dimension") 
    return(W)
  }
  
  W = as.matrix(W)
  D = nrow(W)
  nc = ncol(W)
  isPos = (W>0)
  isNeg = (W<0)
  nPos=colSums(isPos)
  nNeg=colSums(isNeg)
  valPos=sqrt(nNeg/(nPos*(nPos+nNeg)))
  valNeg=-sqrt(nPos/(nNeg*(nPos+nNeg)))
  WW = isPos * outer(rep(1,D),valPos) + 
    isNeg * outer(rep(1,D),valNeg)
  return(WW)
}

#####################################################################################
###############   BENCHMARKING FUNCTIONS  ###########################################
#####################################################################################
# Functions
# 1- pcaRep(x,n) : returns lower dimension(n) of otutable(x) by pca. x must be non-zero.
# 2- pbaRep(x,n) : returns lower dimension(n) of otutable(x) by pba. x must be non-zero.
# 3- amalgamRep(x,n) : returns lower dimension(n) of otutable(x) by amalgamation. x must be non-zero.
# 4- distalBalRep(x,labels,n) : returns lower dimension(n) of otutable(x) by distal balances.
# 5- addLabels(x,labels) : adds a label column to x such as (cirrhosis/healthy) and returns it. 
#                          The resulting dataset is used for Logistic Regression.   
#####################################################################################
## x must be non-zero.
pcaRep <- function(X,numberOfDim) {
  prc <- prcomp(x=log(X),retx=TRUE,rank=numberOfDim,center=TRUE)
  pca.data <- as.data.frame(prc$x) # contains the new principal (ilr) coordinates
  return(pca.data) 
}

#####################################################################################
## x must be non-zero.
pbaRep <- function(x,numberOfDim) {
  library(balance)
  modelPba <- pba(x)
  pb.data <- as.data.frame(modelPba@pba)
  pb.data <- pb.data[,1:(numberOfDim)]
  return(pb.data) 
}

#####################################################################################
## x must be non-zero.
amalgamRep <- function(x,numberOfDim) {
  library(amalgam)
  modelAmalgam <- amalgam(x, n.amalgams = numberOfDim, z = NULL, 
                          objective = objective.keepDist,
                          weight = weight.Nto1)
  amalgam.data<-as.data.frame(modelAmalgam$amalgams)    
  return(amalgam.data) 
}

#####################################################################################
distalBalRep <- function(x,numberOfDim, labels) {
  library(balance)
  sbp <- sbp.fromADBA(x, labels) # get discriminant balances
  sbp <- sbp.subset(sbp) # get distal balances only
  
  modelDistal <- balance.fromSBP(
    x=x, # the data to recast
    y = sbp # the SBP to use
  ) 

  distalBal.data <- as.data.frame(modelDistal[,(ncol(modelDistal)-numberOfDim+1):ncol(modelDistal)])
  
  return(distalBal.data) 
}
#####################################################################################
addLabel<-function(x,labels){
  
  LR.data<-cbind(x,labels)
  names(LR.data)[ncol(LR.data)]<-"label"
  LR.data$label<-as.factor(LR.data$label)
  return(LR.data)
}
