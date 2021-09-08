pdf(file="figS5.pdf", width=7, height=8, paper="a4")

par(mfrow=c(3,2))
par(mar=c(5,5,2,1))
xlab="Position relative to TSS"
ylab="% Occurrence frequency"
ylim=c(0,3.5)
r1=40:183

x=read.table("hg19_c1_yy1_w16.dat")[r1,1]
y=read.table("hg19_c1_yy1_w16.dat")[r1,2]
legend="Class 1"
plot(x,y,type="l",ylim=ylim, xlab=xlab, frame.plot=F,
ylab=ylab, cex.main=1.2, cex.lab=1.3, cex.axis=1.3)
legend("topleft",legend=legend,bty="n", lty=1, cex=1.3)

x=read.table("hg19_c2_yy1_w16.dat")[r1,1]
y=read.table("hg19_c2_yy1_w16.dat")[r1,2]
legend="Class 2"
plot(x,y,type="l",ylim=ylim, xlab=xlab, frame.plot=F,
ylab=ylab, cex.main=1.2, cex.lab=1.3, cex.axis=1.2)
legend("topleft",legend=legend,bty="n", lty=1, cex=1.3)

x=read.table("hg19_c3_yy1_w16.dat")[r1,1]
y=read.table("hg19_c3_yy1_w16.dat")[r1,2]
legend="Class 3"
plot(x,y,type="l",ylim=ylim, xlab=xlab, frame.plot=F,
ylab=ylab, cex.main=1.2, cex.lab=1.3, cex.axis=1.2)
legend("topleft",legend=legend,bty="n", lty=1, cex=1.3)

x=read.table("hg19_c4-6_yy1_w16.dat")[r1,1]
y=read.table("hg19_c4-6_yy1_w16.dat")[r1,2]
legend="Class 4-6"
plot(x,y,type="l",ylim=ylim, xlab=xlab, frame.plot=F,
ylab=ylab, cex.main=1.2, cex.lab=1.3, cex.axis=1.2)
legend("topleft",legend=legend,bty="n", lty=1, cex=1.3)

x=read.table("hg19_tata_yy1_w16.dat")[r1,1]
y=read.table("hg19_tata_yy1_w16.dat")[r1,2]
legend="with TATA"
plot(x,y,type="l",ylim=ylim, xlab=xlab, frame.plot=F,
ylab=ylab, cex.main=1.2, cex.lab=1.3, cex.axis=1.2)
legend("topleft",legend=legend,bty="n", lty=1, cex=1.5)

x=read.table("hg19_notata_yy1_w16.dat")[r1,1]
y=read.table("hg19_notata_yy1_w16.dat")[r1,2]
legend="w/o TATA"
plot(x,y,type="l",ylim=ylim, xlab=xlab, frame.plot=F,
ylab=ylab, cex.main=1.2, cex.lab=1.3, cex.axis=1.2)
legend("topleft",legend=legend,bty="n", lty=1, cex=1.5)

dev.off()

