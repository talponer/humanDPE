dir="./results/"; name="hg19_epd05_C10_s14";

C  =as.matrix(read.table(paste0(dir,name,"_s0.dat"),header=F))

C01=as.matrix(read.table(paste0(dir,name,"_s1.dat"),header=F))
C02=as.matrix(read.table(paste0(dir,name,"_s2.dat"),header=F))
C03=as.matrix(read.table(paste0(dir,name,"_s3.dat"),header=F))
C04=as.matrix(read.table(paste0(dir,name,"_s4.dat"),header=F))
C05=as.matrix(read.table(paste0(dir,name,"_s5.dat"),header=F))
C06=as.matrix(read.table(paste0(dir,name,"_s6.dat"),header=F))
C07=as.matrix(read.table(paste0(dir,name,"_s7.dat"),header=F))
C08=as.matrix(read.table(paste0(dir,name,"_s8.dat"),header=F))
C09=as.matrix(read.table(paste0(dir,name,"_s9.dat"),header=F))
C10=as.matrix(read.table(paste0(dir,name,"_s10.dat"),header=F))

C=C[order(C[1:6,1], decreasing=T),]
bs=matrix(0,nrow=10,ncol=6)
bs[ 1,]=apply(cor(t(C),t(C01)),1,max)
bs[ 2,]=apply(cor(t(C),t(C02)),1,max)
bs[ 3,]=apply(cor(t(C),t(C03)),1,max)
bs[ 4,]=apply(cor(t(C),t(C04)),1,max)
bs[ 5,]=apply(cor(t(C),t(C05)),1,max)
bs[ 6,]=apply(cor(t(C),t(C06)),1,max)
bs[ 7,]=apply(cor(t(C),t(C07)),1,max)
bs[ 8,]=apply(cor(t(C),t(C08)),1,max)
bs[ 9,]=apply(cor(t(C),t(C09)),1,max)
bs[10,]=apply(cor(t(C),t(C10)),1,max) 

colnames(bs)=c(
   sprintf("motif-1 (%4.1f%%)", 100*C[1,1]/sum(C[,1])),
   sprintf("motif-2 (%4.1f%%)", 100*C[2,1]/sum(C[,1])),
   sprintf("motif-3 (%4.1f%%)", 100*C[3,1]/sum(C[,1])),
   sprintf("motif-4 (%4.1f%%)", 100*C[4,1]/sum(C[,1])),
   sprintf("motif-5 (%4.1f%%)", 100*C[5,1]/sum(C[,1])),
   sprintf("motif-6 (%4.1f%%)", 100*C[6,1]/sum(C[,1])))
rownames(bs)=c(
   "bootstrap-01",
   "bootstrap-02",
   "bootstrap-03",
   "bootstrap-04",
   "bootstrap-05",
   "bootstrap-06",
   "bootstrap-07",
   "bootstrap-08",
   "bootstrap-09",
   "bootstrap-10")

library(corrplot)
pdf("./figS2.pdf",8,5)
color <- colorRampPalette(c("red","red","red2","yellow2","green4"))
corrplot(t(bs),method="color", cl.lim=c(0,1),col=color(200), 
   tl.col="black", tl.srt=45, addCoef.col = "black", 
   addgrid.col="white", mar=c(1,1,1,1))
# title(name)
dev.off()

