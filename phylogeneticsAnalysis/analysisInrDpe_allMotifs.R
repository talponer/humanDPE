######################################################################
### This script uses an expectation-maximisation algorithm to separate
### sequences based on their sequence composition.

source("em_functions.R")
## install.packages('seqLogo') or BiocManager::install("seqLogo")
library("seqLogo")
## install.packages('ape')
library('ape') # to do the Neighbour-Joining clustering


# load the data for human for a region from -100 to 100 (to get the
# all CPEs)
hs.file <- "data/Hs_EPDnew_005_hg38.fa"

hs.data <- read.fasta(hs.file)
# Get the INR-box region:
hs.data.clean <- hs.data[-1*unique(which(hs.data == 'N', arr.ind=T)[,1]),]
hs.data.inr <- t(apply(hs.data.clean[,99:132], 1, dna.to.int))

## Set the seeds:
set.seed(1234)

# classify using 3 class (plus the BG)
K <- 3
N <- dim(hs.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
hs.motif <- maximization(hs.data.inr, p, q)

# em
ITER <- 100
for(i in 1:ITER) {
  p <- expectation(hs.data.inr, hs.motif, q)
  q <- colMeans(p)
  hs.motif <- maximization(hs.data.inr, p, q)
  hs.motif <- normalize(hs.motif)
  print(q)
}
# Logo

seqLogo(hs.motif[,,1])
seqLogo(hs.motif[,,2])
seqLogo(hs.motif[,,3])

# amplify the signal:
seqLogo(skew(hs.motif[,,1],3))
seqLogo(skew(hs.motif[,,2],3))
seqLogo(skew(hs.motif[,,3],3))

hs.inr <- hs.motif[,,2]

hs.m1 <- hs.motif[,,order(q, decreasing=T)[1]]
hs.m2 <- hs.motif[,,order(q, decreasing=T)[2]]
hs.m3 <- hs.motif[,,order(q, decreasing=T)[3]]
hs.q <- q

######################################################################
### To make sure that the signal seen here in the DPE region is not
### due to contamination from coding regions, I cleaned the promoter
### selection by pruning all promoters that have a CDS start site in
### the region 0 to +40:
hs.file <- "data/Hs_EPDnew_005_hg38_noCDS.fa"

hs.data <- read.fasta(hs.file)
# Get the INR-box region:
hs.data.clean <- hs.data[-1*unique(which(hs.data == 'N', arr.ind=T)[,1]),]
hs.data.inr <- t(apply(hs.data.clean[,99:132], 1, dna.to.int))

# classify using 3 class (plus the BG)
K <- 3
N <- dim(hs.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
hs.motif <- maximization(hs.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(hs.data.inr, hs.motif, q)
  q <- colMeans(p)
  hs.motif <- maximization(hs.data.inr, p, q)
  hs.motif <- normalize(hs.motif)
  print(q)
}
# Logo

seqLogo(hs.motif[,,1])
seqLogo(hs.motif[,,2])
seqLogo(hs.motif[,,3])

# amplify the signal:
seqLogo(skew(hs.motif[,,1],3))
seqLogo(skew(hs.motif[,,2],3))
seqLogo(skew(hs.motif[,,3],3))

hs.m1 <- hs.motif[,,order(q, decreasing=T)[1]]
hs.m2 <- hs.motif[,,order(q, decreasing=T)[2]]
hs.m3 <- hs.motif[,,order(q, decreasing=T)[3]]
hs.q <- q

## Results are indeed robust and deliver even a better picture of the
## 'DPE' motif

# M. musculus:
mm.file <- "data/Mm_EPDnew_002_mm9.fa"

mm.data <- read.fasta(mm.file)
mm.data.clean <- mm.data[-1*unique(which(mm.data == 'N', arr.ind=T)[,1]),]
mm.data.inr <- t(apply(mm.data.clean[,99:132], 1, dna.to.int))

# classify using 3 class (plus the BG)
K <- 3
N <- dim(mm.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
mm.motif <- maximization(mm.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(mm.data.inr, mm.motif, q)
  q <- colMeans(p)
  mm.motif <- maximization(mm.data.inr, p, q)
  mm.motif <- normalize(mm.motif)
  print(q)
}
# Logo

seqLogo(mm.motif[,,1])
seqLogo(mm.motif[,,2])
seqLogo(mm.motif[,,3])

# amplify the signal:
seqLogo(skew(mm.motif[,,1],2))
seqLogo(skew(mm.motif[,,2],2))
seqLogo(skew(mm.motif[,,3],2))

mm.inr <- mm.motif[,,3]

mm.m1 <- mm.motif[,,order(q, decreasing=T)[1]]
mm.m2 <- mm.motif[,,order(q, decreasing=T)[2]]
mm.m3 <- mm.motif[,,order(q, decreasing=T)[3]]
mm.q <- q


######################################################################
## Now for the promoters with no CDS near the TSS
mm.nocds.file <- "data/Mm_EPDnew_002_mm9_noCDS.fa"

mm.nocds.data <- read.fasta(mm.nocds.file)
# Get the INR-box region:
mm.nocds.data.clean <- mm.nocds.data[-1*unique(which(mm.nocds.data == 'N', arr.ind=T)[,1]),]
mm.nocds.data.inr <- t(apply(mm.nocds.data.clean[,99:132], 1, dna.to.int))

# classify using 3 class (plus the BG)
K <- 3
N <- dim(mm.nocds.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
mm.nocds.motif <- maximization(mm.nocds.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(mm.nocds.data.inr, mm.nocds.motif, q)
  q <- colMeans(p)
  mm.nocds.motif <- maximization(mm.nocds.data.inr, p, q)
  mm.nocds.motif <- normalize(mm.nocds.motif)
  print(q)
}
# Logo

seqLogo(mm.nocds.motif[,,1])
seqLogo(mm.nocds.motif[,,2])
seqLogo(mm.nocds.motif[,,3])

# amplify the signal:
seqLogo(skew(mm.nocds.motif[,,1],3))
seqLogo(skew(mm.nocds.motif[,,2],3))
seqLogo(skew(mm.nocds.motif[,,3],3))

# DPE is still there
mm.m1 <- mm.nocds.motif[,,order(q, decreasing=T)[1]]
mm.m2 <- mm.nocds.motif[,,order(q, decreasing=T)[2]]
mm.m3 <- mm.nocds.motif[,,order(q, decreasing=T)[3]]
mm.q <- q



######################################################################
# D. rerio:
dr.file <- "data/Dr_EPDnew_001_danRer7.fa"

dr.data <- read.fasta(dr.file)
dr.data.clean <- dr.data[-1*unique(which(dr.data == 'N', arr.ind=T)[,1]),]
dr.data.inr <- t(apply(dr.data.clean[,99:132], 1, dna.to.int))

# classify using 3 class (plus the BG)
K <- 3
N <- dim(dr.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
dr.motif <- maximization(dr.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(dr.data.inr, dr.motif, q)
  q <- colMeans(p)
  dr.motif <- maximization(dr.data.inr, p, q)
  dr.motif <- normalize(dr.motif)
  print(q)
}
# Logo

seqLogo(dr.motif[,,1])
seqLogo(dr.motif[,,2])
seqLogo(dr.motif[,,3])

# amplify the signal:
seqLogo(skew(dr.motif[,,1],2))
seqLogo(skew(dr.motif[,,2],2))
seqLogo(skew(dr.motif[,,3],2))

dr.inr <- dr.motif[,,1]

######################################################################
# D. rerio without CDS near TSS:
dr.nocds.file <- "data/Dr_EPDnew_001_danRer7_noCDS.fa"

dr.nocds.data <- read.fasta(dr.nocds.file)
dr.nocds.data.clean <- dr.nocds.data[-1*unique(which(dr.nocds.data == 'N', arr.ind=T)[,1]),]
dr.nocds.data.inr <- t(apply(dr.nocds.data.clean[,99:132], 1, dna.to.int))

# classify using 3 class (plus the BG)
K <- 3
N <- dim(dr.nocds.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
dr.nocds.motif <- maximization(dr.nocds.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(dr.nocds.data.inr, dr.nocds.motif, q)
  q <- colMeans(p)
  dr.nocds.motif <- maximization(dr.nocds.data.inr, p, q)
  dr.nocds.motif <- normalize(dr.nocds.motif)
  print(q)
}
# Logo

seqLogo(dr.nocds.motif[,,1])
seqLogo(dr.nocds.motif[,,2])
seqLogo(dr.nocds.motif[,,3])

# amplify the signal:
seqLogo(skew(dr.nocds.motif[,,1],2))
seqLogo(skew(dr.nocds.motif[,,2],2))
seqLogo(skew(dr.nocds.motif[,,3],2))

dr.nocds.inr <- dr.nocds.motif[,,3]

# In this case the signal is low (but was low even before)
dr.m1 <- dr.nocds.motif[,,order(q, decreasing=T)[1]]
dr.m2 <- dr.nocds.motif[,,order(q, decreasing=T)[2]]
dr.m3 <- dr.nocds.motif[,,order(q, decreasing=T)[3]]
dr.q <- q

######################################################################
# D. melanogaster:
dm.file <- "data/Dm_EPDnew_005_dm6.fa"

dm.data <- read.fasta(dm.file)
#dm.data.clean <- dm.data[-1*unique(which(dm.data == 'N', arr.ind=T)[,1]),]
dm.data.inr <- t(apply(dm.data[,99:132], 1, dna.to.int))

# classify using 3 class (plus the BG)
K <- 3
N <- dim(dm.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
dm.motif <- maximization(dm.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(dm.data.inr, dm.motif, q)
  q <- colMeans(p)
  dm.motif <- maximization(dm.data.inr, p, q)
  dm.motif <- normalize(dm.motif)
  print(q)
}
# Logo

seqLogo(dm.motif[,,1])
seqLogo(dm.motif[,,2])
seqLogo(dm.motif[,,3])

# amplify the signal:
seqLogo(skew(dm.motif[,,1],2))
seqLogo(skew(dm.motif[,,2],2))
seqLogo(skew(dm.motif[,,3],2))

dm.inr <- dm.motif[,,3]

dm.m1 <- dm.motif[,,order(q, decreasing=T)[1]]
dm.m2 <- dm.motif[,,order(q, decreasing=T)[2]]
dm.m3 <- dm.motif[,,order(q, decreasing=T)[3]]
dm.q <- q



######################################################################
# A. mellifera:
am.file <- "data/Am_EPDnew_001_amel5.fa"

am.data <- read.fasta(am.file)
am.data.clean <- am.data[-1*unique(which(am.data == 'N', arr.ind=T)[,1]),]
am.data.inr <- t(apply(am.data.clean[,99:132], 1, dna.to.int))

# classify using 3 class (plus the BG)
K <- 3
N <- dim(am.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
am.motif <- maximization(am.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(am.data.inr, am.motif, q)
  q <- colMeans(p)
  am.motif <- maximization(am.data.inr, p, q)
  am.motif <- normalize(am.motif)
  print(q)
}
# Logo

seqLogo(am.motif[,,1])
seqLogo(am.motif[,,2])
seqLogo(am.motif[,,3])

am.inr <- am.motif[,,2]

am.m1 <- am.motif[,,order(q, decreasing=T)[1]]
am.m2 <- am.motif[,,order(q, decreasing=T)[2]]
am.m3 <- am.motif[,,order(q, decreasing=T)[3]]
am.q <- q


## Interestingly A. mellifera has a DPE-like sequence as well (14% of
## promoters). This is a very strong motif that is slightly different
## from D. melanogaster. There is another large group of promoters
## (42%) that have a canonical Inr but no DPE

######################################################################
# C. elegans:
ce.file <- "data/Ce_EPDnew_001_ce6.fa"

ce.data <- read.fasta(ce.file)
#ce.data.clean <- ce.data[-1*unique(which(ce.data == 'N', arr.ind=T)[,1]),]
ce.data.inr <- t(apply(ce.data[,99:132], 1, dna.to.int))

# classify using 3 class (plus the BG)
K <- 3
N <- dim(ce.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
ce.motif <- maximization(ce.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(ce.data.inr, ce.motif, q)
  q <- colMeans(p)
  ce.motif <- maximization(ce.data.inr, p, q)
  ce.motif <- normalize(ce.motif)
  print(q)
}
# Logo

seqLogo(ce.motif[,,1])
seqLogo(ce.motif[,,2])
seqLogo(ce.motif[,,3])


# amplify the signal:
seqLogo(skew(ce.motif[,,1],2))
seqLogo(skew(ce.motif[,,2],2))
seqLogo(skew(ce.motif[,,3],2))

ce.inr <- ce.motif[,,1]

## C. elegans does not have a DPE or DPE-like element as seen in human
## or mouse (maybe a very weak one)
ce.m1 <- ce.motif[,,order(q, decreasing=T)[1]]
ce.m2 <- ce.motif[,,order(q, decreasing=T)[2]]
ce.m3 <- ce.motif[,,order(q, decreasing=T)[3]]
ce.q <- q


######################################################################
# A. thaliana:
at.file <- "data/At_EPDnew_002_araTha1.fa"

at.data <- read.fasta(at.file)
at.data.clean <- at.data[-1*unique(which(at.data == 'N', arr.ind=T)[,1]),]
at.data.inr <- t(apply(at.data.clean[,99:132], 1, dna.to.int))


# classify using 3 class (plus the BG)
K <- 3
N <- dim(at.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
at.motif <- maximization(at.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(at.data.inr, at.motif, q)
  q <- colMeans(p)
  at.motif <- maximization(at.data.inr, p, q)
  at.motif <- normalize(at.motif)
  print(q)
}
# Logo

seqLogo(at.motif[,,1])
seqLogo(at.motif[,,2])
seqLogo(at.motif[,,3])

at.inr <- at.motif[,,2]

## There is no sign of a DPE element
at.m1 <- at.motif[,,order(q, decreasing=T)[1]]
at.m2 <- at.motif[,,order(q, decreasing=T)[2]]
at.m3 <- at.motif[,,order(q, decreasing=T)[3]]
at.q <- q


######################################################################
# Z. mays:
zm.file <- "data/Zm_EPDnew_001_zm3.fa"

zm.data <- read.fasta(zm.file)
zm.data.clean <- zm.data[-1*unique(which(zm.data == 'N', arr.ind=T)[,1]),]
zm.data.inr <- t(apply(zm.data.clean[,99:132], 1, dna.to.int))


# classify using 3 class (plus the BG)
K <- 3
N <- dim(zm.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
zm.motif <- maximization(zm.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(zm.data.inr, zm.motif, q)
  q <- colMeans(p)
  zm.motif <- maximization(zm.data.inr, p, q)
  zm.motif <- normalize(zm.motif)
  print(q)
}
# Logo

seqLogo(zm.motif[,,1])
seqLogo(zm.motif[,,2])
seqLogo(zm.motif[,,3])

# amplify the signal:
seqLogo(skew(zm.motif[,,1],2))
seqLogo(skew(zm.motif[,,2],2))
seqLogo(skew(zm.motif[,,3],2))

zm.inr <- zm.motif[,,2]

## No sign of a DPE element as well
zm.m1 <- zm.motif[,,order(q, decreasing=T)[1]]
zm.m2 <- zm.motif[,,order(q, decreasing=T)[2]]
zm.m3 <- zm.motif[,,order(q, decreasing=T)[3]]
zm.q <- q

######################################################################
# S. cerevisiae:
sc.file <- "data/Sc_EPDnew_002_sacCer3.fa"

sc.data <- read.fasta(sc.file)
#sc.data.clean <- sc.data[-1*unique(which(sc.data == 'N', arr.ind=T)[,1]),]
sc.data.inr <- t(apply(sc.data[,99:132], 1, dna.to.int))


# classify using 3 class (plus the BG)
K <- 3
N <- dim(sc.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
sc.motif <- maximization(sc.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(sc.data.inr, sc.motif, q)
  q <- colMeans(p)
  sc.motif <- maximization(sc.data.inr, p, q)
  sc.motif <- normalize(sc.motif)
  print(q)
}
# Logo

seqLogo(sc.motif[,,1])
seqLogo(sc.motif[,,2])
seqLogo(sc.motif[,,3])

# amplify the signal:
seqLogo(skew(sc.motif[,,1],2))
seqLogo(skew(sc.motif[,,2],2))
seqLogo(skew(sc.motif[,,3],2))

sc.inr <- sc.motif[,,3]

## No sign of a DPE here eider
sc.m1 <- sc.motif[,,order(q, decreasing=T)[1]]
sc.m2 <- sc.motif[,,order(q, decreasing=T)[2]]
sc.m3 <- sc.motif[,,order(q, decreasing=T)[3]]
sc.q <- q

######################################################################
# S. pombe:
sp.file <- "data/Sp_EPDnew_001_spo2.fa"

sp.data <- read.fasta(sp.file)
#sp.data.clean <- sp.data[-1*unique(which(sp.data == 'N', arr.ind=T)[,1]),]
sp.data.inr <- t(apply(sp.data[,99:132], 1, dna.to.int))


# classify using 3 class (plus the BG)
K <- 3
N <- dim(sp.data.inr)[1]
p <- matrix(0, nrow=N, ncol=K)
# randomly initialise the clusses (using a beta):
for(i in 1:K) {
  p[,i] <- rbeta(N, N**-0.5, 1)
}
p <- p/rowSums(p)
q <- colMeans(p)
sp.motif <- maximization(sp.data.inr, p, q)

# em

# ITER <- 100
for(i in 1:ITER) {
  p <- expectation(sp.data.inr, sp.motif, q)
  q <- colMeans(p)
  sp.motif <- maximization(sp.data.inr, p, q)
  sp.motif <- normalize(sp.motif)
  print(q)
}
# Logo

seqLogo(sp.motif[,,1])
seqLogo(sp.motif[,,2])
seqLogo(sp.motif[,,3])

## # New motif:
## 0.037590411 0.287979090 0.12057858 0.5538519156
## 0.004829764 0.186369714 0.01892506 0.7898754637
## 0.326182815 0.004722191 0.66831023 0.0007847665
## 0.187699784 0.404075992 0.19624289 0.2119813379
## 0.118967301 0.703722036 0.06094556 0.1163650996
## 0.259565682 0.414565927 0.20100145 0.1248669440
## 0.332207274 0.061897559 0.13266831 0.4732268610
## 0.138045171 0.372760107 0.20940541 0.2797893081
## 0.350832805 0.210695350 0.06190165 0.3765701927
## 0.379941861 0.088650314 0.45508435 0.0763234781
## 0.034147994 0.902514835 0.04806313 0.0152740441
## 0.232224856 0.588255671 0.12298792 0.0565315527
## 0.612137364 0.131577618 0.08717611 0.1691089072
## 0.471165459 0.202980132 0.15043762 0.1754167883

## sp.dist <- read.table("data/Sp_newInrDist.dat")

## png("figures/SpNewInrDist.png", width=60, height=60, units="mm", res=600, pointsize=5)
## plot(sp.dist[,1], sp.dist[,2], type="l", lwd=2, frame.plot=F, xlab="Distance from TSS", ylab="Frequency", xlim=c(-50,50))
## points(sp.dist[,1], sp.dist[,3], type="l", lwd=2, col="red")
## dev.off()


# amplify the signal:
seqLogo(skew(sp.motif[,,1],2))
seqLogo(skew(sp.motif[,,2],2))
seqLogo(skew(sp.motif[,,3],2))

sp.inr <- sp.motif[,,3]

## No sign of a DPE but a clear new motif!
sp.m1 <- sp.motif[,,order(q, decreasing=T)[1]]
sp.m2 <- sp.motif[,,order(q, decreasing=T)[2]]
sp.m3 <- sp.motif[,,order(q, decreasing=T)[3]]
sp.q <- q


######################################################################
### Do the NJ clustering with all motifs:
species <- c("Hs","Mm","Dr","Ce","Dm","Am","At","Zm","Sc","Sp")

motif.matrix <- rbind(as.numeric(hs.m1),
                      as.numeric(mm.m1),
                      as.numeric(dr.m1),
                      as.numeric(ce.m1),
                      as.numeric(dm.m1),
                      as.numeric(am.m1),
                      as.numeric(at.m1),
                      as.numeric(zm.m1),
                      as.numeric(sc.m1),
                      as.numeric(sp.m1),
                      as.numeric(hs.m2),
                      as.numeric(mm.m2),
                      as.numeric(dr.m2),
                      as.numeric(ce.m2),
                      as.numeric(dm.m2),
                      as.numeric(am.m2),
                      as.numeric(at.m2),
                      as.numeric(zm.m2),
                      as.numeric(sc.m2),
                      as.numeric(sp.m2),
                      as.numeric(hs.m3),
                      as.numeric(mm.m3),
                      as.numeric(dr.m3),
                      as.numeric(ce.m3),
                      as.numeric(dm.m3),
                      as.numeric(am.m3),
                      as.numeric(at.m3),
                      as.numeric(zm.m3),
                      as.numeric(sc.m3),
                      as.numeric(sp.m3)
                      )

rownames(motif.matrix) <- paste(rep(species, 3),
                                paste('m',
                                      rep(1:3, each=length(species)),
                                      sep=''),
                                sep='.')


motif.dist <- dist(motif.matrix, method="euc")
motif.nj <- nj(motif.dist)

colorsNJ <- c("red","green","cyan","black")
edge.c <- c(rep(1,8), rep(4,1), rep(3,12), rep(4,1), rep(2,14), rep(3,5), rep(1,100))

colorsNJ2 <- c("darkred","darkgreen","blue")
title.c <- c(2,2,1,2,2,2,3,3,3,1,
             1,1,2,1,1,1,1,1,3,1,
             3,3,3,3,2,2,1,3,1,1)

##png("njClusteringAllMotifs_3.png", width=100, height=130, units="mm",
##    res=300, pointsize=6)
plot(motif.nj, type="unrooted", cex=1, main="", font=3, no.margin=T,
     lab4ut="a", edge.col=colorsNJ[edge.c],
     tip.color=colorsNJ2[title.c], label.offset=0.05, rotate.tree=-70)
##dev.off()

pdf("njClusteringAllMotifs_4.pdf", width=10, height=8, pointsize=13)
plot(motif.nj, type="unrooted", cex=1, main="", font=3, no.margin=F,
     lab4ut="a", edge.col=colorsNJ[edge.c],
     tip.color=colorsNJ2[title.c], label.offset=0.05, rotate.tree=85,
     align.tip.label=T)
dev.off()


## Get the consensus for each cluster:

consensus.1 <- (am.m1 + mm.m1 + hs.m1 + dr.m2 + dm.m1 + ce.m1 +
                am.m3 + dm.m3) / 8
consensus.2 <- (sp.m1 + sp.m2 + sc.m3 + sp.m3 + hs.m2 + mm.m2 +
                dr.m1 + dm.m2 + ce.m2 + at.m3 + am.m2 + at.m2 + zm.m2) / 13
consensus.3 <- (at.m1 + zm.m3 + dr.m3 + zm.m1 + sc.m1 + mm.m3 +
                hs.m3 + sc.m2 + ce.m3) / 9

png('consensus1_new.png', width=200, height=100, units="mm", res=200,
     pointsize=6)
seqLogo(skew(consensus.1,2))
dev.off()

png('consensus2_new.png', width=200, height=100, units="mm", res=200,
    pointsize=6)
seqLogo(skew(consensus.2,2))
dev.off()

png('consensus3_new.png', width=200, height=100, units="mm", res=200,
    pointsize=6)
seqLogo(skew(consensus.3,2))
dev.off()

## Make the pdf
pdf('consensus1_new.pdf', width=10, height=5)
seqLogo(skew(consensus.1,2))
dev.off()

pdf('consensus2_new.pdf', width=10, height=5)
seqLogo(skew(consensus.2,2))
dev.off()

pdf('consensus3_new.pdf', width=10, height=5)
seqLogo(skew(consensus.3,2))
dev.off()

######################################################################
### Save the session
######################################################################
## save.image('2020-04-15_promoterAnalysis.RData')
load('2020-04-15_promoterAnalysis.RData')
