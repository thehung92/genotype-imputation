#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
#
INPUT="../data/APMRA96.chr20_QC1.vcf.gz"
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT\n'"
COMMAND=paste(QUERY, INPUT)
df0 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt"))
df0 <- df0 %>%
  mutate(ID=paste(chrom,pos,ref,alt, sep=":"))
#
INPUT="../data/test920.chr20_rmuniq.vcf.gz"
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT\n'"
COMMAND=paste(QUERY, INPUT)
df1<- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt"))
df1 <- df1 %>%
  mutate(ID=paste(chrom,pos,ref,alt, sep=":"))
#
sum(df0$ID %in% df1$ID)
