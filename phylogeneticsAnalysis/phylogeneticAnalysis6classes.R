# partitioning of human promoters with new algorithm

# functions

library('ape')
library('seqLogo')
source("em_functions.R")
source("em_ext_functions.R")

# algorithm parameters for stage 1

seed <- 1      # random number seed
K <- 6         # number of classes
a <- 0.01      # parameter of beta distribution for random seeding
w <- 0.01      # pseudo-weight for base probabilites
s <- 1.4       # motif over-skewing (favors infrequent but highy
               # skewed classes)
x <- 0.02      # pseudo-weight for class probabilites 
nmax <- 300    # total number of motifs to be collected
ITER <- 200    # number iteration per partitioning round
h <- 0.5       # threshold height for cutting tree into subclasses


# stage 1: collection of motifs by multiple rounds of partitioning
#-----------------------------------------------------------------

runEM1 <- function(data, K, a, w, s, x, nmax, iter, seed){
    set.seed(seed)
    Q <- numeric()
    motifs <- list()
    n <- 0
    while(n < nmax) {
        cat(n,"\n")
        ## Initialisation
        N <- dim(data)[1]
        comp <- maximization(data, as.matrix(rep(1, N)), 1)
        p <- matrix(0, nrow=N, ncol=K)
        for(i in 1:K) {p[,i] <- rbeta(N,a,1)}
        p <- p/rowSums(p)
        q <- colMeans(p)
        c <- maximization(data,p,q) 
        c <- normalize(c,w,comp)
        ## em
        for(I in 1:ITER) {
            for(i in 1:K) {c[,,i] <- skew(c[,,i], s)}
            q <- (q+x)/sum(q+x)
            p <- expectation(data, c, q)
            q <- colMeans(p)
            c <- maximization(data, p, q)
            c <- normalize(c, w, comp)
            ## cat(q, "\n")
        }    
        for(i in 1:K) {
            n <- n+1
            motifs[[n]] <- c[,,i]
            Q <- c(Q, q[i])
        }
    }
    ## export data
    L <- length(motifs)
    M <- matrix(nrow=0, ncol=4*dim(data)[2])
    for(i in 1:L) {M <- rbind(M,as.vector(motifs[[i]]))}
    results <- list(M, Q)    
    return(results)
}

runEM2 <- function(motifs, K=6, h, plotTree=FALSE, dpe=NULL){
    M <- motifs[[1]]
    Q <- motifs[[2]]
    tree <- hclust(dist(M))
    if(plotTree){
        if(!is.null(dpe)){
            M1 <- as.data.frame(rbind(M, DPE=as.numeric(dpe)))
            rownames(M1) <- c(paste('m', 1:dim(M)[1], sep=''), 'DpeMotif')
        }else{
            M1 <- as.data.frame(M)
            rownames(M1) <- paste('m', 1:dim(M)[1], sep='')
        }
        treep <- hclust(dist(M1))
        plot(treep)
        abline(h=h, lty=2)
    }
    classes <- cutree(tree, h=h)
    l <- max(classes)
    size <- 1:l
    for(i in 1:l) {size[i] <- length(which(classes == i))}
    ## retain N most frequent classes (one iteration)
    Motifs <- list(); n <- 0;
    for(i in order(size, decreasing=T)[1:K]) {
        n <- n+1;
        x <- matrix(M[which(classes == i),],nrow=size[i]) # necessary
                                                          # if size[i]
                                                          # = 1
        motif <- matrix(colMeans(x),nrow=4)
        Motifs[[n]] <- motif
    }
    P <- 1:K;
    n <- 0
    for(i in order(size, decreasing=T)[1:K]) {
        n <- n+1;
        P[n] <- mean(Q[which(classes == i)])
    }
    P <- P/sum(P)
    ## Re-order classes 
    o <- order(P,decreasing=T)
    P <- P[o]
    M <- Motifs
    for(i in 1:K) {Motifs[[i]] <- M[[o[i]]]}
    ## Output results
    return(M)
}


runEM2.1 <- function(motifs, K=6, k=10, plotTree=FALSE, dpe=NULL){
    M <- motifs[[1]]
    Q <- motifs[[2]]
    tree <- hclust(dist(M))
    if(plotTree){
        if(!is.null(dpe)){
            M1 <- as.data.frame(rbind(M, DPE=as.numeric(dpe)))
            rownames(M1) <- c(paste('m', 1:dim(M)[1], sep=''), 'DpeMotif')
        }else{
            M1 <- as.data.frame(M)
            rownames(M1) <- paste('m', 1:dim(M)[1], sep='')
        }
        treep <- hclust(dist(M1))
        plot(treep)
    }
    classes <- cutree(tree, k=k)
    l <- max(classes)
    size <- 1:l
    for(i in 1:l) {size[i] <- length(which(classes == i))}
    ## retain N most frequent classes (one iteration)
    Motifs <- list(); n <- 0;
    for(i in order(size, decreasing=T)[1:K]) {
        n <- n+1;
        x <- matrix(M[which(classes == i),], nrow=size[i]) ## necessary
                                                           ## if
                                                           ## size[i]
                                                           ## = 1
        motif <- matrix(colMeans(x), nrow=4)
        Motifs[[n]] <- motif
    }
    P <- 1:K;
    n <- 0
    for(i in order(size, decreasing=T)[1:K]) {
        n <- n+1;
        P[n] <- mean(Q[which(classes == i)])
    }
    P <- P/sum(P)
    ## Re-order classes 
    o <- order(P,decreasing=T)
    P <- P[o]
    M <- Motifs
    for(i in 1:K) {Motifs[[i]] <- M[[o[i]]]}
    ## Output results
    return(M)
}


######################################################################
## start with the fly data:

dm.file <- "data/Dm_EPDnew_005_dm6.fa"
dm.data <- read.fasta(dm.file)
# Get the INR-box region:
## dm.data.clean <- dm.data[-1*unique(which(dm.data == 'N', arr.ind=T)[,1]),]
dm.data.inr <- t(apply(dm.data[,91:141], 1, dna.to.int))

dm.motifs1 <- runEM1(dm.data.inr, K, a, w, s, x, nmax, iter, seed)
## dm.FinalMotifs <- runEM2(dm.motifs1, K, h, plotTree=TRUE)
dm.FinalMotifs <- runEM2.1(dm.motifs1, K=K, plotTree=TRUE)
seqLogo(dm.FinalMotifs[[1]])
seqLogo(dm.FinalMotifs[[2]])
seqLogo(dm.FinalMotifs[[3]])
seqLogo(dm.FinalMotifs[[4]])
seqLogo(dm.FinalMotifs[[5]])
seqLogo(dm.FinalMotifs[[6]])



dpeMotif <- dm.FinalMotifs[[6]]

### Analysis of human promoters

fname <- "hg19_epd05.dat"
DATA <- as.matrix(read.table(fname, header=F, row.names=1))
N <- dim(DATA)[1]
L <- dim(DATA)[2]
for(i in 1:N) {for(j in 1:L) {
   if(DATA[i,j] == 0) {DATA[i,j] <- sample(1:4,1) }}}
# select region:
hs.data.inr <- DATA[,41:91]


## hs.file <- "data/Hs_EPDnew_005_hg38.fa"

## hs.data <- read.fasta(hs.file)
## # Get the INR-box region:
## hs.data.clean <- hs.data[-1*unique(which(hs.data == 'N', arr.ind=T)[,1]),]
## hs.data.inr <- t(apply(hs.data.clean[,91:141], 1, dna.to.int))


hs.motifs1 <- runEM1(hs.data.inr, K, a, w, s, x, nmax, iter, seed)
##hs.FinalMotifs <- runEM2(hs.motifs1, K, h, plotTree=TRUE, dpe=dpeMotif)
hs.FinalMotifs <- runEM2.1(hs.motifs1, K=K, plotTree=TRUE, dpe=dpeMotif)
seqLogo(hs.FinalMotifs[[1]])
seqLogo(hs.FinalMotifs[[2]])
seqLogo(hs.FinalMotifs[[3]])
seqLogo(hs.FinalMotifs[[4]])
seqLogo(hs.FinalMotifs[[5]])
seqLogo(hs.FinalMotifs[[6]])



## Mouse data:

mm.file <- "data/Mm_EPDnew_002_mm9.fa"
mm.data <- read.fasta(mm.file)
# Get the INR-box region:
mm.data.clean <- mm.data[-1*unique(which(mm.data == 'N', arr.ind=T)[,1]),]
mm.data.inr <- t(apply(mm.data.clean[,91:141], 1, dna.to.int))

mm.motifs1 <- runEM1(mm.data.inr, K, a, w, s, x, nmax, iter, seed)
##mm.FinalMotifs <- runEM2(mm.motifs1, K=K, h=0.5, plotTree=TRUE, dpe=dpeMotif)
mm.FinalMotifs <- runEM2.1(mm.motifs1, K=K, k=10, plotTree=TRUE, dpe=dpeMotif)
seqLogo(mm.FinalMotifs[[1]])
seqLogo(mm.FinalMotifs[[2]])
seqLogo(mm.FinalMotifs[[3]])
seqLogo(mm.FinalMotifs[[4]])
seqLogo(mm.FinalMotifs[[5]])
seqLogo(mm.FinalMotifs[[6]])


## Zebrafish data:

dr.file <- "data/Dr_EPDnew_001_danRer7.fa"
dr.data <- read.fasta(dr.file)
# Get the INR-box region:
dr.data.clean <- dr.data[-1*unique(which(dr.data == 'N', arr.ind=T)[,1]),]
dr.data.inr <- t(apply(dr.data.clean[,91:141], 1, dna.to.int))

dr.motifs1 <- runEM1(dr.data.inr, K, a, w, s, x, nmax, iter, seed)
## dr.FinalMotifs <- runEM2(dr.motifs1, K=K, h=2.5, plotTree=TRUE, dpe=dpeMotif)
dr.FinalMotifs <- runEM2.1(dr.motifs1, K=K, k=10, plotTree=TRUE, dpe=dpeMotif)
seqLogo(dr.FinalMotifs[[1]])
seqLogo(dr.FinalMotifs[[2]])
seqLogo(dr.FinalMotifs[[3]])
seqLogo(dr.FinalMotifs[[4]])
seqLogo(dr.FinalMotifs[[5]])
seqLogo(dr.FinalMotifs[[6]])


## Bee data:

am.file <- "data/Am_EPDnew_001_amel5.fa"
am.data <- read.fasta(am.file)
# Get the INR-box region:
am.data.clean <- am.data[-1*unique(which(am.data == 'N', arr.ind=T)[,1]),]
am.data.inr <- t(apply(am.data.clean[,91:141], 1, dna.to.int))

am.motifs1 <- runEM1(am.data.inr, K, a, w, s, x, nmax, iter, seed)
## am.FinalMotifs <- runEM2(am.motifs1, K=K, h=2.85, plotTree=TRUE, dpe=dpeMotif)
am.FinalMotifs <- runEM2.1(am.motifs1, K=K, k=10, plotTree=TRUE, dpe=dpeMotif)
seqLogo(am.FinalMotifs[[1]])
seqLogo(am.FinalMotifs[[2]])
seqLogo(am.FinalMotifs[[3]])
seqLogo(am.FinalMotifs[[4]])
seqLogo(am.FinalMotifs[[5]])
seqLogo(am.FinalMotifs[[6]])


## C. elegans data:

ce.file <- "data/Ce_EPDnew_001_ce6.fa"
ce.data <- read.fasta(ce.file)
# Get the INR-box region:
## ce.data.clean <- ce.data[-1*unique(which(ce.data == 'N', arr.ind=T)[,1]),]
ce.data.inr <- t(apply(ce.data[,91:141], 1, dna.to.int))

ce.motifs1 <- runEM1(ce.data.inr, K, a, w, s, x, nmax, iter, seed)
## ce.FinalMotifs <- runEM2(ce.motifs1, K=K, h=1.7, plotTree=TRUE, dpe=dpeMotif)
de.FinalMotifs <- runEM2.1(ce.motifs1, K=K, k=10, plotTree=TRUE, dpe=dpeMotif)
seqLogo(ce.FinalMotifs[[1]])
seqLogo(ce.FinalMotifs[[2]])
seqLogo(ce.FinalMotifs[[3]])
seqLogo(ce.FinalMotifs[[4]])
seqLogo(ce.FinalMotifs[[5]])
seqLogo(ce.FinalMotifs[[6]])


## S. cerevisiae data:

sc.file <- "data/Sc_EPDnew_002_sacCer3.fa"
sc.data <- read.fasta(sc.file)
# Get the INR-box region:
## sc.data.clean <- sc.data[-1*unique(which(sc.data == 'N', arr.ind=T)[,1]),]
sc.data.inr <- t(apply(sc.data[,91:141], 1, dna.to.int))

sc.motifs1 <- runEM1(sc.data.inr, K, a, w, s, x, nmax, iter, seed)
## sc.FinalMotifs <- runEM2(sc.motifs1, K=K, h=3.2, plotTree=TRUE, dpe=dpeMotif)
sc.FinalMotifs <- runEM2.1(sc.motifs1, K=K, k=10, plotTree=TRUE, dpe=dpeMotif)
seqLogo(sc.FinalMotifs[[1]])
seqLogo(sc.FinalMotifs[[2]])
seqLogo(sc.FinalMotifs[[3]])
seqLogo(sc.FinalMotifs[[4]])
seqLogo(sc.FinalMotifs[[5]])
seqLogo(sc.FinalMotifs[[6]])


## S. pombe darta:

sp.file <- "data/Sp_EPDnew_001_spo2.fa"
sp.data <- read.fasta(sp.file)
# Get the INR-box region:
## sp.data.clean <- sp.data[-1*unique(which(sp.data == 'N', arr.ind=T)[,1]),]
sp.data.inr <- t(apply(sp.data[,91:141], 1, dna.to.int))

sp.motifs1 <- runEM1(sp.data.inr, K=K, a, w, s, x, nmax, iter, seed)
## sp.FinalMotifs <- runEM2(sp.motifs1, K=K, h=3.15, plotTree=TRUE, dpe=dpeMotif)
sp.FinalMotifs <- runEM2.1(sp.motifs1, K=K, k=10, plotTree=TRUE, dpe=dpeMotif)
seqLogo(sp.FinalMotifs[[1]])
seqLogo(sp.FinalMotifs[[2]])
seqLogo(sp.FinalMotifs[[3]])
seqLogo(sp.FinalMotifs[[4]])
seqLogo(sp.FinalMotifs[[5]])
seqLogo(sp.FinalMotifs[[6]])


## Z. mays data

zm.file <- "data/Zm_EPDnew_001_zm3.fa"
zm.data <- read.fasta(zm.file)
# Get the INR-box region:
zm.data.clean <- zm.data[-1*unique(which(zm.data == 'N', arr.ind=T)[,1]),]
zm.data.inr <- t(apply(zm.data.clean[,91:141], 1, dna.to.int))

zm.motifs1 <- runEM1(zm.data.inr, K, a, w, s, x, nmax, iter, seed)
## zm.FinalMotifs <- runEM2(zm.motifs1, K=K, h=1.5, plotTree=TRUE, dpe=dpeMotif)
zm.FinalMotifs <- runEM2.1(zm.motifs1, K=K, k=10, plotTree=TRUE, dpe=dpeMotif)
seqLogo(zm.FinalMotifs[[1]])
seqLogo(zm.FinalMotifs[[2]])
seqLogo(zm.FinalMotifs[[3]])
seqLogo(zm.FinalMotifs[[4]])
seqLogo(zm.FinalMotifs[[5]])
seqLogo(zm.FinalMotifs[[6]])

## A. thaliana data

at.file <- "data/At_EPDnew_002_araTha1.fa"
at.data <- read.fasta(at.file)
# Get the INR-box region:
at.data.clean <- at.data[-1*unique(which(at.data == 'N', arr.ind=T)[,1]),]
at.data.inr <- t(apply(at.data.clean[,91:141], 1, dna.to.int))

at.motifs1 <- runEM1(at.data.inr, K, a, w, s, x, nmax, iter, seed)
## at.FinalMotifs <- runEM2(at.motifs1, K=K, h=1.5, plotTree=TRUE, dpe=dpeMotif)
at.FinalMotifs <- runEM2.1(at.motifs1, K=K, k=10, plotTree=TRUE, dpe=dpeMotif)
seqLogo(at.FinalMotifs[[1]])
seqLogo(at.FinalMotifs[[2]])
seqLogo(at.FinalMotifs[[3]])
seqLogo(at.FinalMotifs[[4]])
seqLogo(at.FinalMotifs[[5]])
seqLogo(at.FinalMotifs[[6]])

## 8:14,22:24,27:30,37:43

motifMatrix <- rbind(as.numeric(hs.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(hs.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(hs.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(hs.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(hs.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(hs.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(mm.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(mm.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(mm.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(mm.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(mm.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(mm.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(dr.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(dr.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(dr.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(dr.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(dr.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(dr.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(dm.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(dm.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(dm.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(dm.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(dm.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(dm.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(am.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(am.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(am.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(am.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(am.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(am.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(ce.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(ce.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(ce.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(ce.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(ce.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(ce.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(at.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(at.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(at.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(at.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(at.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(at.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(zm.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(zm.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(zm.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(zm.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(zm.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(zm.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(sc.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(sc.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(sc.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(sc.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(sc.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(sc.FinalMotifs[[6]])[c(8:43)],
                     ##
                     as.numeric(sp.FinalMotifs[[1]])[c(8:43)],
                     as.numeric(sp.FinalMotifs[[2]])[c(8:43)],
                     as.numeric(sp.FinalMotifs[[3]])[c(8:43)],
                     as.numeric(sp.FinalMotifs[[4]])[c(8:43)],
                     as.numeric(sp.FinalMotifs[[5]])[c(8:43)],
                     as.numeric(sp.FinalMotifs[[6]])[c(8:43)]
                     )

species <- c("Hs","Mm","Dr","Dm","Am","Ce","At","Zm","Sc","Sp")
numb <- rep(1:6, length(species))
numb[3] <- 4
numb[4] <- 3

rownames(motifMatrix) <- paste(rep(species, each=6),
                               paste('m', numb, sep=''),
                               sep='.')


motif.dist <- dist(motifMatrix, method="euc")
motif.nj <- nj(motif.dist)

## colorsNJ <- c("red","green","cyan","black")
## edge.c <- c(rep(1,8), rep(4,1), rep(3,12), rep(4,1), rep(2,14), rep(3,5), rep(1,100))

colorsNJ2 <- c("darkred","darkgreen","blue")
title.c <- c(2,2,2,1,3,2,
             2,1,3,3,2,3,
             2,1,1,3,1,1,
             1,2,1,1,1,1,
             2,1,1,1,1,1,
             2,2,3,1,3,1,
             2,2,2,3,2,3,
             2,2,2,2,2,2,
             2,2,2,2,3,2,
             2,3,3,3,2,3)

png("njClusteringAllMotifs_6.png", width=200, height=100, units="mm",
    res=300, pointsize=6)
plot(motif.nj, type="unrooted", cex=1, main="", font=3, no.margin=T,
     lab4ut="a", label.offset=0.05, tip.color=colorsNJ2[title.c], rotate.tree=250)
dev.off()

pdf("njClusteringAllMotifs_6.pdf", width=16, height=8, pointsize=13)
plot(motif.nj, type="unrooted", cex=1, main="", font=3, no.margin=F,
     lab4ut="a", tip.color=colorsNJ2[title.c], label.offset=0.05,
     rotate.tree=250, align.tip.label=T)
dev.off()

######################################################################
## Get the consensus for each group:

allMotifList <- c(hs.FinalMotifs, mm.FinalMotifs, dr.FinalMotifs,
                  dm.FinalMotifs, am.FinalMotifs, ce.FinalMotifs,
                  at.FinalMotifs, zm.FinalMotifs, sc.FinalMotifs,
                  sp.FinalMotifs)

motifCons <- list()
for(I in 1:3){
    mType <- which(title.c == I)
    mAverage <- allMotifList[[mType[1]]]
    for (K in mType[-1]){
        mAverage <- mAverage + allMotifList[[K]]
    }
    mAverage <- mAverage / length(mType)
    motifCons[[I]] <- mAverage
}


pdf("consensus1.pdf", width=12, height=5)
seqLogo(motifCons[[1]])
dev.off()

pdf("consensus2.pdf", width=12, height=5)
seqLogo(motifCons[[2]])
dev.off()

pdf("consensus3.pdf", width=12, height=5)
seqLogo(motifCons[[3]])
dev.off()



######################################################################
### Save image
##save.image('phylogeneticAnalysis6classes.RData')
load('phylogeneticAnalysis6classes.RData')





