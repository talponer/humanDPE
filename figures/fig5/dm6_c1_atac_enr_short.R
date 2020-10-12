infile="dm6_c1_atac_enr.txt"

outfile="dm6_c1_atac_enr_short.pdf"

pdf(file=outfile, width=6, height=4.8)
par(mar=c(5,5,2,1))
data<-read.table(infile, skip=1)[51:91,]
barplot(data[,2], space=c(0,0), col="darkgrey", 
   xlab="Position relative to TSS", xaxt="n",
   ylab="ATAC-seq signal",
   names.arg=data[,1], ylim=c(0,20),
   cex.axis=1.2, cex.lab=1.4, cex.sub=1.5, cex.name=1.4)
points <- length(data[,1])
range <- data[points, 1]+1
step2 <- round(range/10)
axis(1, at=seq(0.5, points, step2),
     labels=seq(data[1,1]+1, data[points, 1]+1, step2), cex.axis=1.2)
legend('topright', bty='n',
   legend=c('D. melanogaster', 'Class 1: Inr+DPE', '37.2%'),
   cex=1.2, adj=c(0, 0.5), col='black')
dev.off()

