######################################################################
### This script uses an expectation-maximisation algorithm to separate
### sequences based on their sequence composition.  Practically it is
### able to find core promoter elements and separate sequences with or
### without them. I will perform the analysis on all EPDnew databases
### restricting the promoters to one per gene (most representative
### promoter).  For the moment the function performs the em on both
### strands (not really needed for promoters) without shifting. I will
### be using Philipp's functions to compare results with Romains'.

source("philipp_em.R")
## install.packages('seqLogo')
library("seqLogo")
## BiocManager::install("seqLogo")
library('seqinr')
## install.packages('ape')
library('ape') # to do the Neighbour-Joining clustering



# load the data for human for a region from -100 to 100 (to get the
# all CPEs)
hs.file <- "data/Hs_EPDnew_005_hg38.fa"

hs.data <- read.fasta(hs.file)
# Get the INR-box region:
hs.data.clean <- hs.data[-1*unique(which(hs.data == 'N', arr.ind=T)[,1]),]
hs.data.inr <- t(apply(hs.data.clean[,99:130], 1, dna.to.int))

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

######################################################################
### To make sure that the signal seen here in the DPE region is not
### due to contamination from coding regions, I cleaned the promoter
### selection by pruning all promoters that have a CDS start site in
### the region 0 to +40:
hs.file <- "data/Hs_EPDnew_005_hg38_noCDS.fa"

hs.data <- read.fasta(hs.file)
# Get the INR-box region:
hs.data.clean <- hs.data[-1*unique(which(hs.data == 'N', arr.ind=T)[,1]),]
hs.data.inr <- t(apply(hs.data.clean[,99:130], 1, dna.to.int))

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

## Results are indeed robust and deliver even a better picture of the
## 'DPE' motif

# M. musculus:
mm.file <- "data/Mm_EPDnew_002_mm9.fa"

mm.data <- read.fasta(mm.file)
mm.data.clean <- mm.data[-1*unique(which(mm.data == 'N', arr.ind=T)[,1]),]
mm.data.inr <- t(apply(mm.data.clean[,99:130], 1, dna.to.int))

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

ITER <- 100
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

######################################################################
## Now for the promoters with no CDS near the TSS
mm.nocds.file <- "data/Mm_EPDnew_002_mm9_noCDS.fa"

mm.nocds.data <- read.fasta(mm.nocds.file)
# Get the INR-box region:
mm.nocds.data.clean <- mm.nocds.data[-1*unique(which(mm.nocds.data == 'N', arr.ind=T)[,1]),]
mm.nocds.data.inr <- t(apply(mm.nocds.data.clean[,99:130], 1, dna.to.int))

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

ITER <- 100
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

######################################################################
# D. rerio:
dr.file <- "data/Dr_EPDnew_001_danRer7.fa"

dr.data <- read.fasta(dr.file)
dr.data.clean <- dr.data[-1*unique(which(dr.data == 'N', arr.ind=T)[,1]),]
dr.data.inr <- t(apply(dr.data.clean[,99:130], 1, dna.to.int))

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

ITER <- 100
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
dr.nocds.data.inr <- t(apply(dr.nocds.data.clean[,99:130], 1, dna.to.int))

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

ITER <- 100
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

######################################################################
# D. melanogaster:
dm.file <- "data/Dm_EPDnew_005_dm6.fa"

dm.data <- read.fasta(dm.file)
#dm.data.clean <- dm.data[-1*unique(which(dm.data == 'N', arr.ind=T)[,1]),]
dm.data.inr <- t(apply(dm.data[,99:130], 1, dna.to.int))

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

ITER <- 100
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


######################################################################
# A. mellifera:
am.file <- "data/Am_EPDnew_001_amel5.fa"

am.data <- read.fasta(am.file)
am.data.clean <- am.data[-1*unique(which(am.data == 'N', arr.ind=T)[,1]),]
am.data.inr <- t(apply(am.data.clean[,99:130], 1, dna.to.int))

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

ITER <- 100
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


## Interestingly A. mellifera has a DPE-like sequence as well (14% of
## promoters). This is a very strong motif that is slightly different
## from D. melanogaster. There is another large group of promoters
## (42%) that have a canonical Inr but no DPE

######################################################################
# C. elegans:
ce.file <- "data/Ce_EPDnew_001_ce6.fa"

ce.data <- read.fasta(ce.file)
#ce.data.clean <- ce.data[-1*unique(which(ce.data == 'N', arr.ind=T)[,1]),]
ce.data.inr <- t(apply(ce.data[,99:130], 1, dna.to.int))

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

ITER <- 100
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

######################################################################
# A. thaliana:
at.file <- "data/At_EPDnew_002_araTha1.fa"

at.data <- read.fasta(at.file)
at.data.clean <- at.data[-1*unique(which(at.data == 'N', arr.ind=T)[,1]),]
at.data.inr <- t(apply(at.data.clean[,99:130], 1, dna.to.int))


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

ITER <- 100
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

######################################################################
# Z. mays:
zm.file <- "data/Zm_EPDnew_001_zm3.fa"

zm.data <- read.fasta(zm.file)
zm.data.clean <- zm.data[-1*unique(which(zm.data == 'N', arr.ind=T)[,1]),]
zm.data.inr <- t(apply(zm.data.clean[,99:130], 1, dna.to.int))


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

ITER <- 100
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

######################################################################
# S. cerevisiae:
sc.file <- "data/Sc_EPDnew_002_sacCer3.fa"

sc.data <- read.fasta(sc.file)
#sc.data.clean <- sc.data[-1*unique(which(sc.data == 'N', arr.ind=T)[,1]),]
sc.data.inr <- t(apply(sc.data[,99:130], 1, dna.to.int))


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

ITER <- 100
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

######################################################################
# S. pombe:
sp.file <- "data/Sp_EPDnew_001_spo2.fa"

sp.data <- read.fasta(sp.file)
#sp.data.clean <- sp.data[-1*unique(which(sp.data == 'N', arr.ind=T)[,1]),]
sp.data.inr <- t(apply(sp.data[,99:130], 1, dna.to.int))


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

ITER <- 100
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

sp.dist <- read.table("data/Sp_newInrDist.dat")

png("figures/SpNewInrDist.png", width=60, height=60, units="mm", res=600, pointsize=5)
plot(sp.dist[,1], sp.dist[,2], type="l", lwd=2, frame.plot=F, xlab="Distance from TSS", ylab="Frequency", xlim=c(-50,50))
points(sp.dist[,1], sp.dist[,3], type="l", lwd=2, col="red")
dev.off()


# amplify the signal:
seqLogo(skew(sp.motif[,,1],2))
seqLogo(skew(sp.motif[,,2],2))
seqLogo(skew(sp.motif[,,3],2))

sp.inr <- sp.motif[,,3]

## No sign of a DPE but a clear new motif!


######################################################################
### Plot all of them with only the INR-box region:

png("figures/HsInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(hs.inr, xfontsize=7, yfontsize=8)
dev.off()
png("figures/HsInrDpeSkew.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(skew(hs.inr,3), xfontsize=7, yfontsize=8)
dev.off()

png("figures/MmInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(mm.inr, xfontsize=7, yfontsize=8)
dev.off()
png("figures/MmInrDpeSkew.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(skew(mm.inr,3), xfontsize=7, yfontsize=8)
dev.off()

png("figures/DrInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(dr.inr, xfontsize=7, yfontsize=8)
dev.off()
png("figures/DrInrDpeSkew.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(skew(dr.inr,3), xfontsize=7, yfontsize=8)
dev.off()

png("figures/DmInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(dm.inr, xfontsize=7, yfontsize=8)
dev.off()

png("figures/AmInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(am.inr, xfontsize=7, yfontsize=8)
dev.off()

png("figures/CeInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(ce.inr, xfontsize=7, yfontsize=8)
dev.off()
png("figures/CeInrDpeSkew.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(ce.inr, xfontsize=7, yfontsize=8)
dev.off()

png("figures/AtInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(at.inr, xfontsize=7, yfontsize=8)
dev.off()

png("figures/ZmInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(zm.inr, xfontsize=7, yfontsize=8)
dev.off()

png("figures/ScInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(sc.inr, xfontsize=7, yfontsize=8)
dev.off()

png("figures/SpInrDpe.png", width=150, height=60, units="mm", res=300, pointsize=8)
seqLogo(sp.inr, xfontsize=7, yfontsize=8)
dev.off()

png("figures/SpNewInr.png", width=100, height=60, units="mm", res=300, pointsize=8)
seqLogo(sp.motif[,1:15,1], xfontsize=7, yfontsize=8)
dev.off()

# Now cluster them (for the momnet exclude A. mellifera and
# S. cerevisiae)

species <- c("H. sapiens","M. musculus","D. rerio","C. elegans","D. melanogaster","A. mellifera","A. thaliana","Z. mays","S. cerevisiae","S. pombe")
#species <- c("H.sapiens","M. musculus","D. rerio","C. elegans","D. melanogaster","A. thaliana","Z. mays","S. pombe")

inr.matrix <- rbind(as.numeric(hs.inr),
                    as.numeric(mm.inr),
                    as.numeric(dr.inr),
                    as.numeric(ce.inr),
                    as.numeric(dm.inr),
                    as.numeric(am.inr),
                    as.numeric(at.inr),
                    as.numeric(zm.inr),
                    as.numeric(sc.inr),
                    as.numeric(sp.inr)
                    )
rownames(inr.matrix) <- species

inr.dist <- dist(inr.matrix, method="euc")
inr.hc <- hclust(inr.dist)

png("figures/inrDpeClustering.png", width=80, height=100, units="mm", res=600, pointsize=7)
plot(inr.hc, main="", xlab="", sub="")
dev.off()

# Now do the Neighbour-Joining clustering:
inr.nj <- nj(inr.dist)


png("figures/inrDpeNjClustering.png", width=100, height=100, units="mm", res=600, pointsize=7)
plot(inr.nj, type="unrooted", cex=1, main="", font=3, no.margin=T, lab4ut="a")
dev.off()


# restrict the clustering to only the Inr region
inr.dist <- dist(inr.matrix[,1:8], method="euc")
inr.hc <- hclust(inr.dist)

png("figures/inrOnlyClustering.png", width=80, height=100, units="mm", res=600, pointsize=7)
plot(inr.hc, main="", xlab="", sub="")
dev.off()

# Now do the Neighbour-Joining clustering:
inr.nj <- nj(inr.dist)


png("figures/inrOnlyNjClustering.png", width=100, height=100, units="mm", res=600, pointsize=7)
plot(inr.nj, type="unrooted", cex=1, main="", font=3, no.margin=T, lab4ut="a")
dev.off()

# restrict the clustering to only the DPE region
inr.dist <- dist(inr.matrix[,9:32], method="euc")
inr.hc <- hclust(inr.dist)

png("figures/dpeOnlyClustering.png", width=80, height=100, units="mm", res=600, pointsize=7)
plot(inr.hc, main="", xlab="", sub="")
dev.off()

# Now do the Neighbour-Joining clustering:
inr.nj <- nj(inr.dist)


png("figures/dpeOnlyNjClustering.png", width=80, height=100, units="mm", res=600, pointsize=7)
plot(inr.nj, type="unrooted", cex=1, main="", font=3, no.margin=T, lab4ut="a")
dev.off()
