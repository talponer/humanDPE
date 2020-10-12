infile="dm6_c3_atac_enr.txt"

outfile="dm6_c3_atac_enr.pdf"

pdf(file=outfile, width=6, height=4.8)
par(mar=c(5,5,2,1))
data<-read.table(infile, skip=1)
barplot(data[,2], space=c(0,0), col="darkgrey", 
   xlab="Position relative to TSS", xaxt="n",
   ylab="ATAC-seq signal",
   names.arg=data[,1], ylim=c(0,20),
   cex.axis=1.2, cex.lab=1.4, cex.sub=1.5, cex.name=1.4)
points <- length(data[,1])
range <- data[points, 1] - data[1,1]
step2 <- round(range/10)
at=0.5+c(0,10,20,30,40,50,60,70,80,90,100)
labels=c(-50,-40,-30,-20,-10,1,11,21,31,41,51)
axis(1, at=at, labels=labels, cex.axis=1.2)
legend('topleft', bty='n',
   legend=c('D. melanogaster', 'Class 3: weak non-canonical Inr', '29.0%'),
   cex=1.2, adj=c(0, 0.5), col='black')
dev.off()

