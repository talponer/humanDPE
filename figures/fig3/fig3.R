library(grid)
source("pwm.R")
source("seqLogo.R")
source("seqLogo2_fig3.R")

pwm=as.matrix(read.table("cons_hg19_K6_s1_1.mat", header=F))
pwm=t(pwm/rowSums(pwm))
pdf("fig3_hg19_class1_logo.pdf",2.7,2.0)
seqLogo2(pwm,xfontsize=6, yfontsize=6, title="Class 1 (weak Inr), 87.9%")
dev.off()

pwm=as.matrix(read.table("cons_hg19_K6_s1_2.mat", header=F))
pwm=t(pwm/rowSums(pwm))
pdf("fig3_hg19_class2_logo.pdf",2.7,2.0)
seqLogo2(pwm,xfontsize=6, yfontsize=6, title="Class 2 (TCT), 4.4%")
dev.off()

pwm=as.matrix(read.table("cons_hg19_K6_s1_3.mat", header=F))
pwm=t(pwm/rowSums(pwm))
pdf("fig3_hg19_class3_logo.pdf",2.7,2.0)
seqLogo2(pwm,xfontsize=6, yfontsize=6, title="Class 3 (Inr + DPE), 3.0%")
dev.off()

pwm=as.matrix(read.table("cons_hg19_K6_s1_4.mat", header=F))
pwm=t(pwm/rowSums(pwm))
pdf("fig3_hg19_class4_logo.pdf",2.7,2.0)
seqLogo2(pwm,xfontsize=6, yfontsize=6, title="Class 4 (Inr + (GCN)n frame 1), 1.9%")
dev.off()

pwm=as.matrix(read.table("cons_hg19_K6_s1_5.mat", header=F))
pwm=t(pwm/rowSums(pwm))
pdf("fig3_hg19_class5_logo.pdf",2.7,2.0)
seqLogo2(pwm,xfontsize=6, yfontsize=6, title="Class 5 (Inr + (GCN)n frame 2), 1.4%")
dev.off()

pwm=as.matrix(read.table("cons_hg19_K6_s1_6.mat", header=F))
pwm=t(pwm/rowSums(pwm))
pdf("fig3_hg19_class6_logo.pdf",2.7,2.0)
seqLogo2(pwm,xfontsize=6, yfontsize=6, title="Class 6 (Inr + (GCN)n frame 3), 1.3%")
dev.off()
