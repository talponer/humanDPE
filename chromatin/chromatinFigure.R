######################################################################
### Script to produce chromatin plots around promoters
library(ggplot2)
library(cowplot)
library(png)
library(grid)

## Human data:
hs.inrdpe.gc <- read.table('../../data/oprof/hg19_epd05_c3_SSSNWWWNSSS.dat')
hs.tata.gc <- read.table('../../data/oprof/hg19_epd05_TATA_SSSNWWWNSSS.dat')
hs.inr.gc <- read.table('../../data/oprof/hg19_epd05_inr_SSSNWWWNSSS.dat')
colnames(hs.inrdpe.gc) <-
    colnames(hs.tata.gc) <-
    colnames(hs.inr.gc) <- c('Position', 'Frequency')

hsGcData <- data.frame(rbind(hs.inrdpe.gc, hs.tata.gc, hs.inr.gc),
                       Type = factor(rep(c('Inr-DPE', 'TATA-box',
                                           'Inr'),
                                         each=dim(hs.inr.gc)[1]),
                                     levels=c('Inr', 'Inr-DPE',
                                             'TATA-box'))
                       )
head(hsGcData)
write.csv(hsGcData, 'humanNucleosomeSignal.csv')

hs1 <- ggplot(hsGcData, aes(x=Position, y=Frequency, col=Type)) +
    geom_line() +
    labs(col='') +
    facet_grid(rows = vars(Type)) +
    annotate("rect", xmin=50, xmax=125, ymin=-Inf, ymax=Inf, alpha=0.1) +
    theme_half_open() +
    panel_border() +
    theme(legend.position = 'none')
ggsave('hs_GCsignal.pdf', hs1, height=7)


## Perform Furier transform
region <- which(hs.inrdpe.gc[,1] > 50 & hs.inrdpe.gc[,1] < 100)
hs.inrdpe.spec <- spec.pgram(hs.inrdpe.gc[region,2], log="yes",
                             plot = FALSE)
hs.tata.spec <- spec.pgram(hs.tata.gc[region,2], log="yes",
                           plot = FALSE)
hs.inr.spec <- spec.pgram(hs.inr.gc[region,2], log="yes",
                          plot = FALSE)
hsSpecData <- data.frame(Length = 1/c(hs.inrdpe.spec$freq,
                                      hs.tata.spec$freq,
                                      hs.inr.spec$freq),
                         Spectrum = c(hs.inrdpe.spec$spec,
                                      hs.tata.spec$spec,
                                      hs.inr.spec$spec),
                         Type = factor(rep(c('Inr-DPE', 'TATA-box',
                                             'Inr'),
                                           each=length(hs.inrdpe.spec$freq)),
                                       levels=c('Inr', 'Inr-DPE',
                                               'TATA-box'))
                         )
head(hsSpecData)
write.csv(hsSpecData, 'humanFourierTransform.csv')


hs2 <- ggplot(hsSpecData, aes(x=Length, y=Spectrum, col=Type)) +
    geom_vline(xintercept=10, col='gray', linetype='dashed') +
    geom_line() +
    labs(col='') +
    scale_x_continuous(breaks=c(0,10,20), limits=c(0,20)) +
    facet_grid(cols = vars(Type)) +
    theme_half_open() +
    panel_border() +
    theme(legend.position = 'none')
ggsave('hs_fourier.pdf', hs2, width=7)
    

## Drosophila data:
dm.inrdpe.gc <- read.table('../../data/oprof/dm6_epd05_c1_SSSNWWWNSSS.dat')
dm.tata.gc <- read.table('../../data/oprof/dm6_epd05_TATA_SSSNWWWNSSS.dat')
dm.inr.gc <- read.table('../../data/oprof/dm6_epd05_c2c3_SSSNWWWNSSS.dat')
colnames(dm.inrdpe.gc) <-
    colnames(dm.tata.gc) <-
    colnames(dm.inr.gc) <- c('Position', 'Frequency')

dmGcData <- data.frame(rbind(dm.inrdpe.gc, dm.tata.gc, dm.inr.gc),
                       Type = factor(rep(c('Inr-DPE', 'TATA-box',
                                           'Inr'),
                                         each=dim(hs.inr.gc)[1]),
                                     levels=c('Inr', 'Inr-DPE',
                                             'TATA-box'))
                       )
head(dmGcData)
write.csv(dmGcData, 'drosophilaNucleosomeSignal.csv')

dm1 <- ggplot(dmGcData, aes(x=Position, y=Frequency, col=Type)) +
    geom_line() +
    labs(col='') +
    facet_grid(rows = vars(Type)) +
    annotate("rect", xmin=50, xmax=125, ymin=-Inf, ymax=Inf, alpha=0.1) +
    theme_half_open() +
    panel_border() +
    theme(legend.position = 'none')
ggsave('dm_GCsignal.pdf', dm1, height=7)

## Perform Furier transform
region <- which(dm.inrdpe.gc[,1] > 50 & dm.inrdpe.gc[,1] < 100)
dm.inrdpe.spec <- spec.pgram(dm.inrdpe.gc[region,2], log="yes",
                             plot = FALSE)
dm.tata.spec <- spec.pgram(dm.tata.gc[region,2], log="yes",
                           plot = FALSE)
dm.inr.spec <- spec.pgram(dm.inr.gc[region,2], log="yes",
                          plot = FALSE)
dmSpecData <- data.frame(Length = 1/c(dm.inrdpe.spec$freq,
                                      dm.tata.spec$freq,
                                      dm.inr.spec$freq),
                         Spectrum = c(dm.inrdpe.spec$spec,
                                      dm.tata.spec$spec,
                                      dm.inr.spec$spec),
                         Type = factor(rep(c('Inr-DPE', 'TATA-box',
                                             'Inr'),
                                           each=length(dm.inrdpe.spec$freq)),
                                       levels=c('Inr', 'Inr-DPE',
                                               'TATA-box'))
                         )
write.csv(dmSpecData, 'drosophilaFourierTransform.csv')

dm2 <- ggplot(dmSpecData, aes(x=Length, y=Spectrum, col=Type)) +
    geom_vline(xintercept=10, col='gray', linetype='dashed') +
    geom_line() +
    xlim(0,20) +
    labs(col='') +
    scale_x_continuous(breaks=c(0,10,20), limits=c(0,20)) +
    facet_grid(cols = vars(Type)) +
    theme_half_open() +
    panel_border() +
    theme(legend.position = 'none')
ggsave('dm_fourier.pdf', dm2, width=7)


######################################################################
### Chromatin data:

## Human nucleosome:
hs.inrdpe.nc <- read.table('../../data/chipcor/nucleosomes/hg19_epd05_c3_vs_Gaffney_all147.dat')
hs.tata.nc <- read.table('../../data/chipcor/nucleosomes/hg19_epd05_TATA_vs_Gaffney_all147.dat')
hs.inr.nc <- read.table('../../data/chipcor/nucleosomes/hg19_epd05_c1_vs_Gaffney_all147.dat')
colnames(hs.inrdpe.nc) <-
    colnames(hs.tata.nc) <-
    colnames(hs.inr.nc) <- c('Position', 'Enrichment')
hsNcData <- data.frame(rbind(hs.inrdpe.nc, hs.tata.nc, hs.inr.nc),
                       Type = factor(rep(c('Inr-DPE', 'TATA-box',
                                           'Inr'),
                                         each=dim(hs.inr.nc)[1]),
                                     levels=c('Inr', 'Inr-DPE',
                                             'TATA-box'))
                       )
head(hsNcData)
write.csv(hsNcData, 'humanMNase.csv')

hs3 <- ggplot(hsNcData, aes(x=Position, y=Enrichment, col=Type)) +
    geom_line() +
    labs(col='') +
    facet_grid(rows = vars(Type)) +
    theme_half_open() +
    panel_border() +
    theme(legend.position = 'none')
ggsave('hs_nucleosomeSignal.pdf', hs3, height=7)

## CAGE
hs.inrdpe.cage <- read.table('../../data/chipcor/CAGE/hg19_epd05_c3_vs_FANTOMall.dat')
hs.tata.cage <- read.table('../../data/chipcor/CAGE/hg19_epd05_TATA_vs_FANTOMall.dat')
hs.inr.cage <- read.table('../../data/chipcor/CAGE/hg19_epd05_c1_vs_FANTOMall.dat')
colnames(hs.inrdpe.cage) <-
    colnames(hs.tata.cage) <-
    colnames(hs.inr.cage) <- c('Position', 'Enrichment')
hsCAGEData <- data.frame(rbind(hs.inrdpe.cage, hs.tata.cage, hs.inr.cage),
                       Type = factor(rep(c('Inr-DPE', 'TATA-box',
                                           'Inr'),
                                         each=dim(hs.inr.cage)[1]),
                                     levels=c('Inr', 'Inr-DPE',
                                             'TATA-box'))
                       )
head(hsCAGEData)
write.csv(hsCAGEData, 'humanCAGE.csv')

hs4 <- ggplot(hsCAGEData, aes(x=Position, y=Enrichment, col=Type)) +
    geom_line() +
    labs(col='') +
    xlim(-50,50) +
    theme_half_open() +
    theme(legend.position=c(0.1, 0.9))
ggsave('hs_CAGEsignal.pdf', hs4)

######################################################################
## Drosophila nucleosome:
dm.inrdpe.nc <- read.table('../../data/chipcor/nucleosomes/dm6_epd05_c1.SRR2038264.dat')
dm.tata.nc <- read.table('../../data/chipcor/nucleosomes/dm6_epd05_TATA.SRR2038264.dat')
dm.inr.nc <- read.table('../../data/chipcor/nucleosomes/dm6_epd05_c2c3.SRR2038264.dat')
colnames(dm.inrdpe.nc) <-
    colnames(dm.tata.nc) <-
    colnames(dm.inr.nc) <- c('Position', 'Enrichment')

dmNcData <- data.frame(rbind(dm.inrdpe.nc, dm.tata.nc, dm.inr.nc),
                       Type = factor(rep(c('Inr-DPE', 'TATA-box',
                                           'Inr'),
                                         each=dim(dm.inr.nc)[1]),
                                     levels=c('Inr', 'Inr-DPE',
                                             'TATA-box'))
                       )
head(dmNcData)
write.csv(dmNcData, 'drosophilaMNase.csv')

dm3 <- ggplot(dmNcData, aes(x=Position, y=Enrichment, col=Type)) +
    geom_line() +
    labs(col='') +
    facet_grid(rows = vars(Type)) +
    theme_half_open() +
    panel_border() +
    theme(legend.position = 'none')
ggsave('dm_nucleosomeSignal.pdf', dm3, height=7)

## CAGE
dm.inrdpe.cage <- read.table('../../data/chipcor/CAGE/dm6_epd05_c1_vs_Machibase.dat')
dm.tata.cage <- read.table('../../data/chipcor/CAGE/dm6_epd05_TATA_vs_Machibase.dat')
dm.inr.cage <- read.table('../../data/chipcor/CAGE/dm6_epd05_c2c3_vs_Machibase.dat')
colnames(dm.inrdpe.cage) <-
    colnames(dm.tata.cage) <-
    colnames(dm.inr.cage) <- c('Position', 'Enrichment')

dmCAGEData <- data.frame(rbind(dm.inrdpe.cage, dm.tata.cage, dm.inr.cage),
                       Type = factor(rep(c('Inr-DPE', 'TATA-box',
                                           'Inr'),
                                         each=dim(dm.inr.cage)[1]),
                                     levels=c('Inr', 'Inr-DPE',
                                             'TATA-box'))
                       )
head(dmCAGEData)
write.csv(dmCAGEData, 'drosophilaCAGE.csv')

dm4 <- ggplot(dmCAGEData, aes(x=Position, y=Enrichment, col=Type)) +
    geom_line() +
    labs(col='') +
    xlim(-50,50) +
    theme_half_open() +
    theme(legend.position=c(0.1, 0.9))
ggsave('dm_CAGEsignal.pdf', dm4)


######################################################################
### Plot all together

pdf('chromatinStates.pdf', height=8, width=17)
plot_grid(hs1, hs2, hs3, hs4, dm1, dm2, dm3, dm4, labels="AUTO",
          ncol=4, label_size=18)
dev.off()

pdf('chromatinStatesAlt.pdf', height=8, width=17)
plot_grid(hs4, hs3, hs1, hs2, dm4, dm3, dm1, dm2, labels="AUTO",
          ncol=4, label_size=18)
dev.off()



######################################################################
### Plot DPE distribution around human and fly promoters:
humanDPE <- read.table("../../data/oprof/hg19_epd05_DPEelement.dat")
flyDPE <- read.table("../../data/oprof/dm6_epd05_DPEelement.dat")
colnames(humanDPE) <- colnames(flyDPE) <- c("Position", "DPE",
                                            "Shuffle")
DPEimg <- readPNG('DPEelement.png')
g <- rasterGrob(DPEimg, interpolate=TRUE)

dpeData <- data.frame(Position = rep(humanDPE[,1], 4),
                     Frequency = c(humanDPE[,2], humanDPE[,3],
                                    flyDPE[,2], flyDPE[,3]),
                      Type = rep(rep(c("DPE", "Shuffle"),
                                     each=dim(humanDPE)[1]), 2),
                      Organism = rep(c("H. sapiens",
                                       "D. melanogaster"),
                                     each=dim(humanDPE)[1]*2)
                      )
head(dpeData)
write.csv(dpeData, 'dpeDistribution.csv')

dp1 <- ggplot(dpeData, aes(x=Position, y=Frequency, col=Type)) +
    geom_line(size=1.3) +
    scale_color_manual(values=c("black", "gray")) +
    labs(col='') +
    facet_grid(rows = vars(Organism)) +
    ## annotate("rect", xmin=50, xmax=125, ymin=-Inf, ymax=Inf,
    ##          alpha=0.1) +
    annotation_custom(g, xmin=-100, xmax=-20, ymin=20, ymax=30) +
    theme_half_open() +
    panel_border()
ggsave('DPE_distribution.pdf', dp1, height=7)
