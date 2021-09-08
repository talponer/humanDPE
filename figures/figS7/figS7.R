skew = function(mat,f) {
  x = mat**f
  mat = t(t(x)/colSums(x))
  return(mat)
}

library(seqLogo)

tollo=as.matrix(read.table("tollo.mat"))
cg10479=as.matrix(read.table("cg10479.mat"))
dm6_c1=as.matrix(read.table("cons_dm6_K3_s1_1.mat")[25:39,])
hg19_c3=as.matrix(read.table("cons_hg19_K6_s1_3.mat")[25:39,])

pdf(file="figS7_tollo_overskewed_by2.pdf"  , width=6, height=4); seqLogo(skew(t(tollo  ),2)); dev.off()
pdf(file="figS7_cg10479_overskewed_by2.pdf", width=6, height=4); seqLogo(skew(t(cg10479),2)); dev.off()
pdf(file="figS7_dm6_c1_overskewed_by2.pdf" , width=6, height=4); seqLogo(skew(t(dm6_c1 ),2)); dev.off()
pdf(file="figS7_hg19_c3_overskewed_by2.pdf", width=6, height=4); seqLogo(skew(t(hg19_c3),2)); dev.off()

# correlation plot 

motifs=as.vector(tollo)
motifs=cbind(motifs,as.vector(cg10479))
motifs=cbind(motifs,as.vector(dm6_c1))
motifs=cbind(motifs,as.vector(hg19_c3))
M=cor(motifs,motifs)
rownames(M)=c("tollo","cg10479", "dm6_c1", "hg19_c3")
colnames(M)=c("tollo","cg10479", "dm6_c1", "hg19_c3")

library(corrplot)
pdf(file="figS7_logo_correlation_plot.pdf")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF",
    "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200), type="upper", diag=T,
    tl.col="black", tl.srt=45, addCoef.col = "black",
    addgrid.col="white", mar=c(6,6,6,6))
dev.off()

