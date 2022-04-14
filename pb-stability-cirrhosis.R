# Studying stability of Principal balances (R-mode stability)
# cirrhosis data after zero substitution (Asli Data)

# read data
Data = read.csv(file="cirrhosis_otuno0.csv")
colnames(Data)
#edit(Data)
Data<-otu.no0
######################################
# computes the Aitchison distance between
# two compositions x,y
Adistxy<-function(x,y){
  xclr = log(x)-mean(log(x))
  yclr = log(y)-mean(log(y))
  Adistxy = sqrt(sum((xclr -yclr)^2))
  return(Adistxy)
}
#####################################

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
  hdist = sqrt(mvariation(x,optnor="none"))
  hdist=as.dist(hdist)
  # cluster
  hmerge=hclust(d=hdist,method="ward.D2")
  # signs from $merge; contrast matrix
  Vsigns=mmerge2sign(hmerge$merge)
  V = mbuildcontrast(Vsigns)
  rownames(V)=colnames(x)
  return(V)
}
#################################################




n = nrow(Data)
d = ncol(Data)-1
n
d
#number of (bootstrap) resamplings
nsam = 50

# sizes of nsizes subsamples
nsubs = c(40,50,60,70,80,90,100,110,120,130)
nsizes=length(nsubs)
nsizes

#original analysis
xor = Data    # original data
Vxor = mPBclustvar(xor)  # r-mode custer
xorilr = milr(xor,Vxor)  # ilr coordinates following Vxor
vxorilr = diag(var(xorilr))  # composition of variances

# init variances of coordinates and distances
# option 0: to neutral element
# option 1: to original ;
optref=0
# init matrix of coordinates
vxilr = array(0,dim=c(nsam,nsizes,d))
# tree distances to original analysis (optref=1), 
# to neutral (optref=0)
Neutral = rep(1,d+1)

# storage of distances 
Tdist = array(0,dim=c(nsam,nsizes,(d+1))) 

# Loop on bootstrap samples
for(isam in 1:nsam){
    # Loop on size of samples
    for(isub in 1:nsizes){ 
      # bootstraping subsample of size ns
      ns = nsubs[isub]
      sam = ceiling(runif(n=ns,min=0,max=n))
      x=Data[sam,]
# Rmode cluster: (n,n-1) contrast matrix
      Vx = mPBclustvar(x)
#      sign(Vx)
# ilr coordinates resp Vx
      xilr = milr(x,Vx)
# var-covar of coordinates
      vxilr[isam,isub,1:d] = diag(var(xilr))
# tree distances for each number of groups
# d is the number of ilr-coordinates equals number of groups-1
      for(igr in 1:d){
       #igr is the number of groups
      if(optref==0){
        Tdist[isam,isub,igr] = Adistxy(Neutral[(1:igr)],
                                    vxilr[isam,isub,(1:igr)])
        }
      if(optref==1){
        Tdist[isam,isub,igr] = Adistxy(vxorilr[1:igr],
                                vxilr[isam,isub,(1:igr)])
        }
      } # igr
    }  # isub
}  # isam


#......................................
# plot number of groups against quantiles of distances
# and original data estimate

# number of groups to plot
numgr = 100
# loop on number of groups
# size of samples is set to the original n=nsubs[nsizes]=130
Q = matrix(0,nrow=numgr,ncol=8)
for(igr in 2:numgr){
  zdist = Tdist[1:nsam,nsizes,igr]
  # quantiles for zdist
  Q[igr,1:7]=quantile(zdist,probs=c(1,0.9,0.75,0.5,0.25,0.1,0))
}
# norms of composition of variances in original sample
dNor = rep(0,numgr)
for(jgroup in 2:numgr){
  dNor[jgroup] = Adistxy(Neutral[1:(jgroup-1)],vxorilr[1:(jgroup-1)])
}

pdf(file="Qnorm-nPMG.pdf",width=6,height=5)
plot(3:numgr,Q[2:(numgr-1),1],col="black",
     xlim=c(0,numgr),ylim=c(0,9),pch=1,cex=0.7,
     xlab="number of PMG's",ylab="quantiles of norms")
points(3:numgr,Q[2:(numgr-1),7],col="black",pch=1,cex=0.7)
lines(3:numgr,dNor[3:100],col="red") 
lines(3:numgr,Q[2:(numgr-1),4],col="blue")
lines(3:numgr,Q[2:(numgr-1),3],col="black")
lines(3:numgr,Q[2:(numgr-1),5],col="black")
dev.off()

# plot of norms/num groups
pdf(file="Qnorm2-nG-nPMG.pdf",width=6,height=5)
plot(3:numgr,(Q[2:(numgr-1),1])^2/3:numgr,col="black",
     xlim=c(0,numgr),ylim=c(0,1),pch=1,cex=0.7,
     xlab="number of PMG's",ylab="quantiles of norms^2/nPMGs")
points(3:numgr,(Q[2:(numgr-1),7])^2/3:numgr,col="black",pch=1,cex=0.7)
lines(3:numgr,(dNor[3:100])^2/3:numgr,col="red") 
lines(3:numgr,(Q[2:(numgr-1),4])^2/3:numgr,col="blue")
lines(3:numgr,(Q[2:(numgr-1),3])^2/3:numgr,col="black")
lines(3:numgr,(Q[2:(numgr-1),5])^2/3:numgr,col="black")
dev.off()


#......................................
# plot quantiles of norms against size of bootstraped samples 
# 
nsubs
nsizes
# reference number of groups
nref =27

QQ = matrix(0,nrow=nsizes,ncol=8)
for(isub in 1:nsizes){
  # sample of norms
  zdist = Tdist[1:isam,isub,nref]
  # compute quantiles
  QQ[isub,1:7] = quantile(zdist,probs=c(1,0.9,0.75,0.5,0.25,0.1,0))
}
# original
dNref= Adistxy(Neutral[1:(nref-1)],vxorilr[1:(nref-1)])
dNref

pdf(file="Qnorms-ssize.pdf")
plot(nsubs[1:nsizes],QQ[1:nsizes,1],col="black",ylim=c(3.5,4.7),
     ylab="quantiles of norm",xlab="size of sample")
points(nsubs[1:nsizes],QQ[1:nsizes,7],col="black")
lines(nsubs[1:nsizes],QQ[1:nsizes,3],col="black")
lines(nsubs[1:nsizes],QQ[1:nsizes,4],col="blue")
lines(nsubs[1:nsizes],QQ[1:nsizes,5],col="black")
lines(nsubs[c(1,nsizes)],c(dNref,dNref),col="red")
dev.off()
