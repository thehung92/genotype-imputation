#!/usr/bin/env Rscript
args=commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<=1) {
  stop("At least two argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==2) {
  # default output dir
  args[3] = "./"
  print(paste("take", args[1:2], "and output in", args[3]))
} else if (length(args)==3) {
  print(paste("take", args[1:2], "and output in", args[3]))
}
# test
# args <- c("/Users/hung/Data/DGV4VN/Impute_protocol/run_4/data/test1014.chr20_filter.vcf.gz", "/Users/hung/Data/DGV4VN/APMRA/analysis/output/APMRA96_annotate.vcf.gz", "./")
#
library(data.table)
library(tidyverse)
#
# read PID from sequencing data
INPUT=args[1]
# INPUT="/Users/hung/Data/DGV4VN/Impute_protocol/run_4/data/test1014.chr20_filter.vcf.gz"
QUERY="source ~/.bashrc; bcftools query -l"
COMMAND=paste(QUERY, INPUT)
df0 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names="PID")
df0 <- df0 %>%
  mutate(PID2=str_extract(PID, "\\d{4}"))
#
# read PID from genotyping data
INPUT=args[2]
# INPUT="/Users/hung/Data/DGV4VN/APMRA/analysis/output/APMRA96_annotate.vcf.gz"
QUERY="source ~/.bashrc; bcftools query -l"
COMMAND=paste(QUERY, INPUT)
df1 <- fread(cmd=COMMAND, header=FALSE,
             col.names="PID")
df1$PID2 <- sapply(df1$PID, function(x){
  gsub(".CEL", "", x) -> x
  sprintf("%04d", as.numeric(x))})
#
OUTDIR=args[3]
# OUTDIR="./"
# write samples list used in VN ref panel
df0 %>%
  filter(! PID2 %in% df1$PID2) %>%
  select(PID) -> output0
OUTPUT0=paste0(OUTDIR, "list-ref-samples.txt")
fwrite(output0, file=OUTPUT0, quote=FALSE, sep="\t",
       row.names = FALSE, col.names = FALSE)
# write samples list used in test input and write bcftools renaming file
df1 %>%
  filter(PID2 %in% df0$PID2) %>%
  select(PID) -> output1
OUTPUT1=paste0(OUTDIR, "list-input-samples.txt")
fwrite(output1, file=OUTPUT1, quote=FALSE, sep="\t",
       row.names = FALSE, col.names = FALSE)
# bcftools renaming file : old name \t new name
OUTPUT2=paste0(OUTDIR, "samples-rename.txt")
df1 %>%
  filter(PID2 %in% df0$PID2) %>%
  left_join(df0, by="PID2") %>%
  select(1,3) -> output2
fwrite(output2, file=OUTPUT2, quote=FALSE, sep="\t",
       row.names = FALSE, col.names = FALSE)
#
print("Done")