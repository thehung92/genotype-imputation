#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
library(ggpubr)
#
#
source("function.test_accuracy.R")
source("function.plot_accuracy.R")
source("function.plot_hqvariants.R")
source("function.class_maf.R")
#
#### import dosage ####
# PID
INPUT="../output/APMRA94.chr20_impute.1KGP3.dose.vcf.gz"
QUERY1="source ~/.bashrc ; bcftools query -l"
COMMAND1=paste(QUERY1, INPUT) # colnames
df1 <- fread(cmd=COMMAND1, header=FALSE, col.names="PID")
# 1kgp
INPUT="../output/APMRA94.chr20_impute.1KGP3.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
ds.1kgp <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
             col.names = c("ID", df1$PID))
# test920
INPUT="../output/APMRA94.chr20_impute.test920.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
ds.vn920 <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
             col.names = c("ID", df1$PID))
# merge
INPUT="../output/APMRA94.chr20_impute.merge-1KGP3-test920.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
ds.merge <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
             col.names = c("ID", df1$PID))
# merge-2
INPUT="../output/APMRA94.chr20_impute.merge-reciprocal.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
ds.merge2 <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
                  col.names = c("ID", df1$PID))
# SG10K
INPUT="../output/APMRA94.chr20_impute.SG10K.dose.vcf.gz"
QUERY2="source ~/.bashrc ; bcftools query -i 'INFO/IMPUTED=1' -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n'"
COMMAND2=paste(QUERY2, INPUT) # dosage
ds.sg10k <- fread(cmd=COMMAND2, sep="\t", header=FALSE,
                  col.names = c("ID", df1$PID))
# true
INPUT="../data/test1014.chr20_phase.vcf.gz"
fwrite(df1, file="samples_list.txt", quote=FALSE, col.names=FALSE)
QUERY="source ~/.bashrc ; bcftools query -S samples_list.txt -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n'"
COMMAND=paste(QUERY, INPUT)
gt.vn94 <- fread(cmd=COMMAND, sep="\t", header=FALSE,
                 col.names = c("ID", df1$PID))
# convert matrix
mat <- as.matrix(gt.vn94)
mat[mat == "0|0"] <- "0"
mat[mat == "0|1" | mat == "1|0"] <- "1"
mat[mat == "1|1"] <- "2"
mat <- as_tibble(mat) %>%
  mutate(across(matches("VN"), as.numeric))

#### import MAF from 1KGP3 ####
MAF_1KGP3 <- fread(file="~/Data/1KGP_LC/frequency/MAF_1KGP3.chr20.txt", sep="\t", header=TRUE)
MAF_1KGP3 %>%
  # cbind(cut(MAF_1KGP3$MAF,c(1e-4,1e-3,5e-3,1e-2,5e-2,0.2,0.5))) %>%
  cbind(cut(MAF_1KGP3$MAF,c(0,1e-3,1e-2,0.1,0.5))) %>%
  as_tibble() -> MAF_1KGP3
colnames(MAF_1KGP3)[3] <- "MAF_group"

#### intersected sites between datasets ####
# create identical intersected sites
# 286698
Reduce(intersect,list(ds.1kgp$ID, ds.vn920$ID, ds.merge$ID, ds.sg10k$ID, gt.vn94$ID)) -> m

#### test accuracy ####
r2_vn920 <- test_accuracy(mat, ds.vn920, MAF_1KGP3)
r2_1kgp <- test_accuracy(mat, ds.1kgp, MAF_1KGP3)
r2_merge <- test_accuracy(mat, ds.merge, MAF_1KGP3)
r2_merge2 <- test_accuracy(mat, ds.merge2, MAF_1KGP3)
r2_sg10k <- test_accuracy(mat, ds.sg10k, MAF_1KGP3)

#### plot accuracy ####
p3 <- plot_accuracy(x=ls(pattern="r2_"), output="accuracy_5-ref-panel.png")
# ggsave(filename="accuracy0.png", plot=p1, units="in", width=6, height=8, dpi=300)



#### import rsq ####
# 1KGP3
rsq_1kgp <- fread(file="../output/APMRA94.chr20_impute.1KGP3.info", sep="\t", header = TRUE) %>%
  select(c(1,7,8)) %>% mutate(across(Genotyped, as.factor))
# test920
rsq_test920 <- fread(file="../output/APMRA94.chr20_impute.test920.info", sep="\t", header = TRUE) %>%
  select(c(1,7,8)) %>% mutate(across(Genotyped, as.factor))
# merge
rsq_merge <- fread(file="../output/APMRA94.chr20_impute.merge-1KGP3-test920.info", sep="\t", header = TRUE) %>%
  select(c(1,7,8)) %>% mutate(across(Genotyped, as.factor))
# merge2
rsq_merge2 <- fread(file="../output/APMRA94.chr20_impute.merge-reciprocal.info", sep="\t", header = TRUE) %>%
  select(c(1,7,8)) %>% mutate(across(Genotyped, as.factor))
# SG10K
rsq_sg10k <- fread(file="../output/APMRA94.chr20_impute.SG10K.info", sep="\t") %>%
  select(c(1,7,8)) %>% mutate(across(Genotyped, as.factor))
#### plot high quality imputed variants ####
p4 <- plot_hqvariants(x=ls(pattern="rsq_"), output="high-quality-variants_5-ref-panel.png")

#
save(list=c("p3", "p4"), file="branch-hgdp-p34.Rdata")
