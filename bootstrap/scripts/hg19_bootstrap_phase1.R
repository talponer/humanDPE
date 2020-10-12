# Functions

source("./scripts/Philipp_em2_functions.R")

# Filenames:

infile="./data/hg19_epd05_s0.dat"
infiles=c(
   "./data/hg19_epd05_s1.dat",
   "./data/hg19_epd05_s2.dat",
   "./data/hg19_epd05_s3.dat",
   "./data/hg19_epd05_s4.dat",
   "./data/hg19_epd05_s5.dat",
   "./data/hg19_epd05_s6.dat",
   "./data/hg19_epd05_s7.dat",
   "./data/hg19_epd05_s8.dat",
   "./data/hg19_epd05_s9.dat",
   "./data/hg19_epd05_s10.dat")

outfile="./results/hg19_epd05_K6_s14_s0.dat"
outfiles=c(
   "./results/hg19_epd05_K6_s14_s1.dat",
   "./results/hg19_epd05_K6_s14_s2.dat",
   "./results/hg19_epd05_K6_s14_s3.dat",
   "./results/hg19_epd05_K6_s14_s4.dat",
   "./results/hg19_epd05_K6_s14_s5.dat",
   "./results/hg19_epd05_K6_s14_s6.dat",
   "./results/hg19_epd05_K6_s14_s7.dat",
   "./results/hg19_epd05_K6_s14_s8.dat",
   "./results/hg19_epd05_K6_s14_s9.dat",
   "./results/hg19_epd05_K6_s14_s10.dat") 

# algorithm parameters for stage 1

seed=1      # random number seed
K=6         # number of classes
a=0.01      # parameter of beta distribution for random seeding
w=0.01      # pseudo-weight for base probabilites
s=1.4       # motif over-skewing (favors infrequent but highy skewed classes)
x=0.02      # pseudo-weight for class probabilites
nmax=300    # total number of motifs to be collected
ITER=200    # number iteration per partitioning round

# stage 1: collection of motifs by multiple rounds of partitioning
#-----------------------------------------------------------------

j=as.integer(commandArgs(trailingOnly=TRUE)[1])

   if(j == 0) {
      data=as.matrix(read.table(infile, header=F, row.names=1))
   } else {
      data=as.matrix(read.table(infiles[j], header=F, row.names=1))}

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

   L=length(motifs)
   M=matrix(nrow=0, ncol=4*dim(data)[2])
   for(i in 1:L) {M=rbind(M,as.vector(motifs[[i]]))}

   if(j == 0) {
      write.table(cbind(Q,M),outfile,
      quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
   } else {
      write.table(cbind(Q,M),outfiles[j],
      quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")}

