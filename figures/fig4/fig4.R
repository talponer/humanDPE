# Fig 4A (fig4_dm6_atac_seq_long.pdf)

dm6_c1=read.table("dm6_c1_atac_enr.txt", skip=1)[1:100,2]; dm6_c1=dm6_c1/mean(dm6_c1)
dm6_c2=read.table("dm6_c2_atac_enr.txt", skip=1)[1:100,2]; dm6_c2=dm6_c2/mean(dm6_c2)
dm6_c3=read.table("dm6_c3_atac_enr.txt", skip=1)[1:100,2]; dm6_c3=dm6_c3/mean(dm6_c3)
x=c(-50:-1,1:50)

outfile="fig4_dm6_atac_seq_long.pdf"
pdf(file=outfile, width=8, height=5)

par(mar=c(5,5,2,1))

col1 <- "#E41A1C" #red
col2 <- "#377EB8" #blue
col3 <- "#4DAF4A" #green

plot(dm6_c1, frame.plot=F, col=col1,
   xlab="Position relative to TSS", xaxt="n",
   ylab="ATAC-seq signal",
   t="l", ylim=c(0,4), lwd=2,
   cex.axis=1.2, cex.lab=1.4, cex.sub=1.5)
points <- c(1,11,21,31,41,51,60,70,80,90,100)
axis(1, at=points, labels=x[points],
   cex.axis=1.2)
lines(dm6_c2, col=col2,  t="l", ylim=c(0,4), lwd=2)
lines(dm6_c3, col=col3, t="l", ylim=c(0,4), lwd=2)
legend("topright", bty="n", lwd=2,
   legend=c("Class 1 (Inr+DPE)","Class 2 (Inr)","Class 3 (weak Inr)"),
   cex=1.2, adj=c(0, 0.5), col=c(col1,col2,col3))

dev.off()

# Fig 4B (fig4_hg19_atac_seq_long.pdf)

hg19_c1=read.table("hg19_c1_atac_enr.txt", skip=1)[1:100,2]; hg19_c1=hg19_c1/mean(hg19_c1)
hg19_c2=read.table("hg19_c2_atac_enr.txt", skip=1)[1:100,2]; hg19_c2=hg19_c2/mean(hg19_c2)
hg19_c3=read.table("hg19_c3_atac_enr.txt", skip=1)[1:100,2]; hg19_c3=hg19_c3/mean(hg19_c3)
x=c(-50:-1,1:50)

outfile="fig4_hg19_atac_seq_long.pdf"
pdf(file=outfile, width=8, height=5)

par(mar=c(5,5,2,1))

col1 <- "#E41A1C" #red
col2 <- "#377EB8" #blue
col3 <- "#4DAF4A" #green

plot(hg19_c1, frame.plot=F, col=col3,
   xlab="Position relative to TSS", xaxt="n",
   ylab="ATAC-seq signal",
   t="l", ylim=c(0,4), lwd=2,
   cex.axis=1.2, cex.lab=1.4, cex.sub=1.5)
points <- c(1,11,21,31,41,51,60,70,80,90,100)
axis(1, at=points, labels=x[points],
   cex.axis=1.2)
lines(hg19_c2, col=col2,  t="l", ylim=c(0,4), lwd=2)
lines(hg19_c3, col=col1, t="l", ylim=c(0,4), lwd=2)
legend("topright", bty="n", lwd=2,
   legend=c("Class 1 (weak Inr)","Class 2 (TCT)","Class 3 (Inr+DPE)"),
   cex=1.2, adj=c(0, 0.5), col=c(col3,col2, col1))

dev.off()

# Fig 4C (fig4_atac_seq_short_dpe.pdf)

dm6_c1=read.table("dm6_c1_atac_enr.txt", skip=1)[51:91,2]
hg19_c3=read.table("hg19_c3_atac_enr.txt", skip=1)[51:91,2]
dm6_c1=dm6_c1/mean(dm6_c1)
hg19_c3=hg19_c3/mean(hg19_c3)
x=1:41

outfile="fig4_atac_seq_short_dpe.pdf"
pdf(file=outfile, width=7, height=5)

col1="grey"; col2="black";

par(mar=c(5,5,2,1))

plot(dm6_c1, frame.plot=F, col=col1,
   xlab="Position relative to TSS", xaxt="n",
   ylab="ATAC-seq signal",
   t="l", ylim=c(0,4), lwd=2,
   cex.axis=1.2, cex.lab=1.4, cex.sub=1.5)
points <- length(x)
range <- x[points] - x[1]
step2 <- round(range/10)
axis(1, at=seq(0.5, points, step2/1), labels=seq(x[1], x[points], step2),
   cex.axis=1.2)
lines(hg19_c3, col=col2, t="l", ylim=c(0,4), lwd=2)
legend('topright', bty="n", lwd=2,
   legend=c("Drosophila Class 1 (Inr+DPE)", "Human Class 3 (Inr+DPE)"),
   cex=1.2, adj=c(0, 0.5), col=c(col1,col2))

dev.off()

# Fig 4D (fig4_atac_correlation_matrix.pdf)

r=51:91
dm6_c1= read.table("dm6_c1_atac_enr.txt")[r,2]
dm6_c2= read.table("dm6_c2_atac_enr.txt")[r,2]
dm6_c3= read.table("dm6_c3_atac_enr.txt")[r,2]
hg19_c1=read.table("hg19_c1_atac_enr.txt")[r,2]
hg19_c2=read.table("hg19_c2_atac_enr.txt")[r,2]
hg19_c3=read.table("hg19_c3_atac_enr.txt")[r,2]

all=cbind(dm6_c1,dm6_c2,dm6_c3,hg19_c1,hg19_c2,hg19_c3)
M=cor(all,all)

library(corrplot)
pdf(file="fig4_atac_correlation_matrix.pdf")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF",
    "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200), type="upper", diag=T,
    tl.col="black", tl.srt=45, addCoef.col = "black",
    addgrid.col="white", mar=c(1,1,1,1))
dev.off()

