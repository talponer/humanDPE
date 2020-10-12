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
pdf(file="atac_correlation_plot.pdf")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF",
    "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200), type="upper", diag=T,
    tl.col="black", tl.srt=45, addCoef.col = "black",
    addgrid.col="white", mar=c(1,1,1,1))
dev.off()

