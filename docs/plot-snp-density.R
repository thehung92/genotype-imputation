#!/usr/bin/env Rscript
# plot snp density
# BiocManager::install("karyoploteR")
#
library(data.table)
library(tidyverse)
library(karyoploteR)
# read data 1KGP
INPUT="../data/1KGP3.chr20_rmuniq.vcf.gz"
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC\n'"
COMMAND=paste(QUERY, INPUT)
df0 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt", "ac"))
grange0 <- toGRanges(df0[,c("chrom", "pos", "pos")])
# read data VN
INPUT="../data/test920.chr20_rmuniq.vcf.gz"
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC\n'"
COMMAND=paste(QUERY, INPUT)
df1 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt", "ac"))
grange1 <- toGRanges(df1[,c("chrom", "pos", "pos")])
# read data SG10K
INPUT="../data/SG10K.chr20_QC3.vcf.gz"
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC\n'"
COMMAND=paste(QUERY, INPUT)
df2 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt", "ac"))
grange2 <- toGRanges(df2[,c("chrom", "pos", "pos")])

# plot chr20 of hg38
OUTDIR="/Users/hung/Data/DGV4VN/Impute_protocol/run_4/docs/"
OUTPUT=paste0(OUTDIR, "temp3.pdf")
pdf(file=OUTPUT, bg="white",
    width=8, height=4)
kp <- plotKaryotype(genome="hg38", chromosomes = "chr20", plot.type=2)
#
kp <- kpPlotDensity(kp, data=grange0, window.size = 100000, r0=0, r1=0.5)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5)
kp <- kpPlotDensity(kp, data=grange2, window.size = 100000, r0=0.5, r1=1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.5, r1=1)
#
kp <- kpPlotDensity(kp, data=grange1, data.panel=2, window.size = 100000)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=2)
dev.off()

# exclude the centromere region based on one interval
INPUT="../data/1KGP3.chr20_rmuniq.vcf.gz"
ARG="/Users/hung/Data/DGV4VN/Impute_protocol/run_1/centromere/centromere.hg38.chr20.bed "
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC\n' -T ^"
COMMAND=paste0(QUERY, ARG, INPUT)
df0 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt", "ac"))
grange0 <- toGRanges(df0[,c("chrom", "pos", "pos")])
#
INPUT="../data/test920.chr20_rmuniq.vcf.gz"
ARG="/Users/hung/Data/DGV4VN/Impute_protocol/run_1/centromere/centromere.hg38.chr20.bed "
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC\n' -T ^"
COMMAND=paste0(QUERY, ARG, INPUT)
df1 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt", "ac"))
grange1 <- toGRanges(df1[,c("chrom", "pos", "pos")])
#
OUTDIR="/Users/hung/Data/DGV4VN/Impute_protocol/run_4/docs/"
OUTPUT=paste0(OUTDIR, "temp2.pdf")
pdf(file=OUTPUT, bg="white",
    width=8, height=4)
kp <- plotKaryotype(genome="hg38", chromosomes = "chr20", plot.type=2)
#
kp <- kpPlotDensity(kp, data=grange0, window.size = 100000)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density)
#
kp <- kpPlotDensity(kp, data=grange1, data.panel=2, window.size = 100000)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=2)
#
dev.off()
#
#### test ####
set.seed(1000)

data <- createRandomRegions(nregions=20000)

kp <- plotKaryotype("hg19", plot.type=2, chromosomes="chr1")

kp <- kpPlotDensity(kp, data)
kpAxis(kp, ymin = 0, ymax=kp$latest.plot$computed.values$max.density)

kp <- kpPlotDensity(kp, data, data.panel=2, col="#CCCCFF",  ymax=20, lwd=2)
kpAxis(kp, ymin = 0, ymax=20, data.panel=2)

kp <- kpLines(kp, data=kp$latest.plot$computed.values$windows, y=kp$latest.plot$computed.values$density, col="black", r0=0.5, r1=1, data.panel=2, ymax=20)



rm(list=ls())
