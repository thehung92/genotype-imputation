#!/usr/bin/env Rscript
args=commandArgs(trailingOnly = TRUE)
#
library(data.table)
library(tidyverse)
#
INPUT=args[1]
# INPUT="../data/HGDP.b38_QC1.vcf.gz"
QUERY="source ~/.bashrc ; bcftools query -f '%CHROM %POS %REF %ALT %ID\n'"
COMMAND=paste(QUERY, INPUT)
df0 <- fread(cmd=COMMAND, header=FALSE, sep=" ",
             col.names=c("chrom", "pos", "ref", "alt", "id"))
#
set.seed(123456)
output <- sample_frac(df0[,1:2], 0.1)
#
OUTDIR=args[2]
OUTPUT=paste0(OUTDIR, "masked.snps.list")
fwrite(output, file=OUTPUT, quote=FALSE, sep="\t", col.names = FALSE)