#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
source("function.convert_gt_matrix.R")
#
#### import dosage ####
# PID
INPUT="../output/HGDP-vn94.chr20_QC1.vcf.gz"
QUERY="source ~/.bashrc ; bcftools query -l"
CMD=paste(QUERY, INPUT) # colnames
vt0 <- fread(cmd=CMD, header=FALSE, col.names="PID") %>% unlist()
# read all imputed file in a list of data frame
registerDoParallel(8)
#
vt1 <-list.files("../output/", pattern="HGDP.*gz") %>% paste0("../output/", .)
ldf <- foreach (i=vt1) %dopar% {
  INPUT=i
  QUERY="source ~/.bashrc ; bcftools query -T ../output/masked.snps.list -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n'"
  CMD=paste(QUERY, INPUT)
  y <- fread(cmd=CMD, header=FALSE, sep="\t",
             col.names=c("ID", vt0))
  return(y)
}
# create intersection of ID from all df
ldf2 <- lapply(ldf, function(x) {
  x[,"ID"] %>% unlist()
})
Reduce(intersect, ldf2) -> vt2

# convert matrix to 0-1-2 genotype, expecting NA value in true genotype matrix
ldf2 <- lapply(ldf, function(x){
  # x is data frame from ldf
  # subset to intersection of variant ID
  # remove ID column
  y <- filter(x, ID %in% vt2) %>%
    select(-1) %>% as.matrix()
  # fill z with condition from y
  y[y=="0/0" | y=="0|0"] <- "0"
  y[y=="0/1" | y=="0|1"] <- "1"
  y[y=="1/0" | y=="1|0"] <- "1"
  y[y=="1/1" | y=="1|1"] <- "2"
  # convert value to numeric
  z <- as_tibble(y) %>%
    mutate_all(as.numeric) %>%
    as.matrix()
  #
  return(z)
})
# rename list element
vt3 <- gsub(".*impute\\.(.+)\\.dose.*", "\\1", x=vt1)
vt3[length(vt3)] <- "true"
names(ldf2) <- vt3
# calculate error rate
true <- ldf2[[6]]
ldf2[[6]] <- NULL
sapply(ldf2, function(x) {
  # x is matrix from ldf2
  sum(!x==true, na.rm=TRUE)
})
#### split ldf2 into subset of population ####
# read info from excel file
INPUT="/Users/hung/Data/HGDP/hgdp_stanford/HGDPid_populations.xls"
df0 <- readxl::read_excel(INPUT, range = "A1:F1065")
# filter asian population
df1 <- df0 %>%
  filter(Id %in% vt0) %>% # subset with available sample in vcf
  filter(Region=="Asia") %>% # subset asian 
  select(1,3) %>%
  mutate_at(2, as.factor)
#### parallel main function ####
registerDoParallel(8)
foreach (i=levels(df1$population), .combine=rbind) %dopar% {
  # subset column based on population in df1
  print(i)
  # get sample ID from df1
  samples <- df1 %>%
    filter(population==i) %>%
    select(1) %>%
    unlist()
  # return(c(i, length(samples)))
  # subset from ldf2 and true
  mat0 <- true[, samples]
  lmat <- lapply(ldf2, function(x) {
    # x is matrix element of ldf2 list
    x[, samples]
  })
  # compute error rate
  output <- sapply(lmat, function(x) {
    # x is matrix from ldf2
    sum(!x==mat0, na.rm=TRUE)
  })
  #
  return(c(i, output))
} -> df.output

#### save Rdata ####
save(list=ls(), file="branch-hgdp.Rdata")
