#!/usr/bin/env Rscript
# check centromere region
#
library(data.table)
library(tidyverse)
# VN
INPUT="/Users/hung/Data/DGV4VN/Impute_protocol/run_4/data/test920.chr20_rmuniq.vcf.gz"
ARG="/Users/hung/Data/DGV4VN/Impute_protocol/run_1/centromere/centromere.hg38.chr20.bed"
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC\n' -T"
COMMAND=paste(QUERY, ARG, INPUT)
df0 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt", "ac"))
#

# 1KGP
INPUT="/Users/hung/Data/DGV4VN/Impute_protocol/run_4/data/1KGP3.chr20_rmuniq.vcf.gz"
ARG="/Users/hung/Data/DGV4VN/Impute_protocol/run_1/centromere/centromere.hg38.chr20.bed"
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT %INFO/AC\n' -T"
COMMAND=paste(QUERY, ARG, INPUT)
df1 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt", "ac"))
#



