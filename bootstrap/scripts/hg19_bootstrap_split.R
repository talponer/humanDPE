# Functions

source("./scripts/Philipp_em2_functions.R")

# sequence input

fname="./data/hg19_epd05.dat"

DATA=as.matrix(read.table(fname, header=F, row.names=1))
N=dim(DATA)[1]; L=dim(DATA)[2]
for(i in 1:N) {for(j in 1:L) {
   if(DATA[i,j] == 0) {DATA[i,j]=sample(1:4,1) }}}

# select region:

data=DATA[,41:91]

   outfile="./data/hg19_epd05_s0.dat"
   write.table(data, outfile, col.names=FALSE, quote=F)

N=dim(data)[1]

J=10;
for(j in 1:J) {
   subset=sample(1:N,N,replace=T)
   subset_data=data[subset,]
   rownames=rep("-",N)
   for(i in 1:N) {rownames[i]=paste0(rownames(subset_data)[i], sprintf("_%i", i))}
   rownames(subset_data)=rownames
   outfile=paste0("data/hg19_epd05_s",sprintf("%i", j),".dat")
   write.table(subset_data, outfile, col.names=FALSE, quote=F)
   }

