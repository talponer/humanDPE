# Functions

source("./scripts/Philipp_em2_functions.R")

# major loop 

for(j in 0:10) {

# Motif input

infile="results/hg19_epd05_K6_s14_s0.dat"
infiles=c(
   "results/hg19_epd05_K6_s14_s1.dat",
   "results/hg19_epd05_K6_s14_s2.dat",
   "results/hg19_epd05_K6_s14_s3.dat",
   "results/hg19_epd05_K6_s14_s4.dat",
   "results/hg19_epd05_K6_s14_s5.dat",
   "results/hg19_epd05_K6_s14_s6.dat",
   "results/hg19_epd05_K6_s14_s7.dat",
   "results/hg19_epd05_K6_s14_s8.dat",
   "results/hg19_epd05_K6_s14_s9.dat",
   "results/hg19_epd05_K6_s14_s10.dat")

if(j == 0) {
   M=read.table(infile)      
} else {
   M=read.table(infiles[j])}

L=dim(M)[1];Q=M[,1];M=M[,2:(dim(M)[2])]
motifs=list(); for(i in 1:L) {
   motifs[[i]]=matrix(as.numeric(M[i,]),nrow=4)}

# algorithm parameters for stage 2

K=10        # number of classes
h=0.5       # height parameter for tree-cut

# clustering

h=0.5
classes=cutree(hclust(dist(M)),h=h)
if(max(classes) <= K) {classes=cutree(hclust(dist(M)),k=K)}
l=max(classes); size=1:l
for(i in 1:l) {size[i] = length(which(classes == i))} 

# retain K most frequent classes (one iteration)

Motifs=list(); n=0;
for(i in order(size, decreasing=T)[1:K]) {n=n+1;
   if(size[i] == 1) {
      Motifs[[n]]=matrix(as.vector(M[which(classes == i),]),nrow=4)
   } else {Motifs[[n]]=matrix(colMeans(M[which(classes == i),]),nrow=4)}}
P=1:K; n=0
for(i in order(size, decreasing=T)[1:K]) {n=n+1;
   P[n]=mean(Q[which(classes == i)])
   }
P=P/sum(P)

# re-order motifs

# o=order(P, decreasing=T)
# motifs=Motifs; for(i in 1:K) {Motifs[[i]]=motifs[[o[i]]]}
# P=P[o]

# write consensus motifs

outfile="results/hg19_epd05_C10_s14_s0.dat"
outfiles=c(
   "results/hg19_epd05_C10_s14_s1.dat",
   "results/hg19_epd05_C10_s14_s2.dat",
   "results/hg19_epd05_C10_s14_s3.dat",
   "results/hg19_epd05_C10_s14_s4.dat",
   "results/hg19_epd05_C10_s14_s5.dat",
   "results/hg19_epd05_C10_s14_s6.dat",
   "results/hg19_epd05_C10_s14_s7.dat",
   "results/hg19_epd05_C10_s14_s8.dat",
   "results/hg19_epd05_C10_s14_s9.dat",
   "results/hg19_epd05_C10_s14_s10.dat")

C=matrix(0,nrow=length(Motifs), ncol=4*dim(Motifs[[1]])[2])
for(i in 1:K) {C[i,]=as.vector(Motifs[[i]])}
if(j == 0) { 
   write.table(cbind(P,C),outfile,
   quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
} else {
   write.table(cbind(P,C),outfiles[j],
   quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")}

# end of major loop

}
