# partitioning of human promoters with new algorithm

# functions

library(seqLogo)
source("Philipp_em2_functions.R")

# Data input

fname="hg19_epd05.dat"

DATA=as.matrix(read.table(fname, header=F, row.names=1))
N=dim(DATA)[1]; L=dim(DATA)[2]
for(i in 1:N) {for(j in 1:L) {
   if(DATA[i,j] == 0) {DATA[i,j]=sample(1:4,1) }}}

# select region:

data=DATA[,41:91]
 
# algorithm parameters for stage 1

seed=1      # random number seed
K=3         # number of classes
a=0.01      # parameter of beta distribution for random seeding
w=0.01      # pseudo-weight for base probabilites
s=1.0       # motif over-skewing (favors infrequent but highy skewed classes)
x=0.02      # pseudo-weight for class probabilites 
nmax=3      # total number of motifs to be collected
ITER=200    # number iteration per partitioning round

name="hg19_K3_s1" # Run-name for output filenames

# stage 1: collection of motifs by multiple rounds of partitioning 
#-----------------------------------------------------------------

set.seed(seed); Q=numeric(); motifs=list(); n=0; 
while(n < nmax) {print(n)

  N=dim(data)[1]
  comp=maximization(data,as.matrix(rep(1,N)),1)
  p=matrix(0,nrow=N, ncol=K)
  for(i in 1:K) {p[,i] = rbeta(N,a,1)}
  p=p/rowSums(p)
  q = colMeans(p)
  c = maximization(data,p,q) 
  c=normalize(c,w,comp)

# em

  for(I in 1:ITER) {
    for(i in 1:K) {c[,,i]=skew(c[,,i],s)}
    q=(q+x)/sum(q+x)  
    p=expectation(data,c,q)
    q=colMeans(p)
    c=maximization(data,p,q)
    c=normalize(c,w,comp)
    print(q)
    }
  for(i in 1:K) {n=n+1; motifs[[n]]=c[,,i]; Q=c(Q,q[i])}
  }

# export data

L=length(motifs)
M=matrix(nrow=0, ncol=4*dim(data)[2])
for(i in 1:L) {M=rbind(M,as.vector(motifs[[i]]))}
write.table(cbind(Q,M),paste0("results/",name,".dat"), 
quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")


# import data (to start from here with saved stage-1 motifs) 

if(F) { # disabled by default 
M=read.table(paste0("results/",name,".dat"))
L=dim(M)[1];Q=M[,1];M=M[,2:(dim(M)[2])]
motifs=list(); for(i in 1:L) {
   motifs[[i]]=matrix(as.numeric(M[i,]),nrow=4)}}

# make images (for visualization of stage-1 motif)

if(F) { # disabled by default 
L=length(motifs)
for(i in 1:L) { 
   motif=paste0(name,sprintf("_%03i_%4.2f", i, Q[i]*100))
   file=paste0("img/",motif,".jpg")
   jpeg(file, width = 540, height = 360, pointsize=12) 
   seqLogo(as.matrix(motifs[[i]])); dev.off()}}

# Stage 2: extract and merge repeatedly found motifs 
# (useless if K=nmax - kept for compatibility reasons)
#-----------------------------------------------------------------------

# Parameters

# K=3; s=1.0 # define again if starting from here
h=0.5   # threshold height for cutting tree into subclasses

# motif clustering 

classes=cutree(hclust(dist(M)),h=h)
l=max(classes); size=1:l
for(i in 1:l) {size[i] = length(which(classes == i))}

# retain N most frequent classes (one iteration)

Motifs=list(); n=0;
for(i in order(size, decreasing=T)[1:K]) {n=n+1;
   x=matrix(M[which(classes == i),],nrow=size[i]) # necessary if size[i] = 1
   motif=matrix(colMeans(x),nrow=4)
   Motifs[[n]]=motif}
P=1:K; n=0
for(i in order(size, decreasing=T)[1:K]) {n=n+1;
   P[n]=mean(Q[which(classes == i)])
   }
P=P/sum(P)

# Re-order classes 

o=order(P,decreasing=T)
P=P[o]
M=Motifs; for(i in 1:K) {Motifs[[i]] = M[[o[i]]]}

# export data (one table by motif)

for(i in 1:K) {M=matrix(sprintf("%6.4f", t(Motifs[[i]])), ncol=4)
write.table(M,paste0("results/cons_",name,sprintf("_%i",i),".mat"),
quote=FALSE, row.names=FALSE, col.names=FALSE)}

# all in one

L=length(Motifs)
MC=matrix(nrow=0, ncol=4*dim(data)[2])
for(i in 1:L) {MC=rbind(MC,as.vector(Motifs[[i]]))}
write.table(cbind(P,MC),paste0("results/cons_",name,".dat"),
quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

# make images

L=length(Motifs)
for(i in 1:L) {  
   motif=paste0("cons_",name,sprintf("_%03i_%4.2f", i, P[i]*100))
   file=paste0("img/",motif,".jpg")
   jpeg(file, width = 540, height = 360, pointsize=12) 
   seqLogo(Motifs[[i]]); dev.off()}

# Assign promoters to classes:

epd=read.table("hg19_epd05.sga")
C = array(data=0, dim=c(4, dim(data)[2],K))
for(i in 1:K) {C[,,i]=skew(Motifs[[i]],s)}
p=expectation(data,C,P)
class=rep(0,dim(p)[1])
for(i in 1:dim(p)[1]) {class[i]=which.max(p[i,])}
write.table(cbind(epd,class),"hg19_epd05_classes.sga",
sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

