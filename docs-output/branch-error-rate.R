#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
library(ggpubr)
#
source("function.convert_gt_matrix.R")
source("function.plot_error_rate.R")
#### import genotype ####
# PID
INPUT="../output/APMRA94.chr20_impute.1KGP3.dose.vcf.gz"
QUERY1="source ~/.bashrc ; bcftools query -l"
COMMAND1=paste(QUERY1, INPUT) # colnames
df1 <- fread(cmd=COMMAND1, header=FALSE, col.names="PID")
# 1kgp
INPUT="../output/APMRA94.chr20_impute.1KGP3.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
gt.1kgp <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
                 col.names = c("ID", df1$PID))
gt.1kgp <- convert_gt_matrix(gt.1kgp)
# test920
INPUT="../output/APMRA94.chr20_impute.test920.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
gt.vn920 <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
                  col.names = c("ID", df1$PID))
gt.vn920 <- convert_gt_matrix(gt.vn920)
# merge
INPUT="../output/APMRA94.chr20_impute.merge-1KGP3-test920.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
gt.merge <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
                  col.names = c("ID", df1$PID))
gt.merge <- convert_gt_matrix(gt.merge)
# merge-2
INPUT="../output/APMRA94.chr20_impute.merge-reciprocal.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
gt.merge2 <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
                   col.names = c("ID", df1$PID))
gt.merge2 <- convert_gt_matrix(gt.merge2)
# SG10K
INPUT="../output/APMRA94.chr20_impute.SG10K.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
gt.sg10k <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
                  col.names = c("ID", df1$PID))
gt.sg10k <- convert_gt_matrix(gt.sg10k)
# true
INPUT="../data/test1014.chr20_phase.vcf.gz"
fwrite(df1, file="samples_list.txt", quote=FALSE, col.names=FALSE)
QUERY="source ~/.bashrc ; bcftools query -S samples_list.txt -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n'"
COMMAND=paste(QUERY, INPUT)
gt.vn94 <- fread(cmd=COMMAND, sep="\t", header=FALSE,
                 col.names = c("ID", df1$PID))
true <- convert_gt_matrix(gt.vn94)
#### import MAF ####
MAF_1KGP3 <- fread(file="~/Data/1KGP_LC/frequency/MAF_1KGP3.chr20.txt", sep="\t", header=TRUE)
MAF_1KGP3 %>%
  cbind(cut(MAF_1KGP3$MAF,c(0,1e-3,1e-2,0.1,0.5))) %>%
  as_tibble() -> MAF_1KGP3
colnames(MAF_1KGP3)[3] <- "MAF_group"
#
df0 <- plot_error_rate(true, gt.1kgp, MAF_1KGP3)
