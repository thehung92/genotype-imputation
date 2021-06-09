#!/usr/bin/env Rscript
args=commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
  print(paste("take", args[1], "and produce", args[2]))
} else if (length(args)==2) {
  print(paste("take", args[1], "and produce", args[2]))
}
#
library(data.table)
library(tidyverse)
#
INPUT=args[1]
# INPUT="/Users/hung/Data/DGV4VN/APMRA/analysis/output/APMRA96_annotate.vcf.gz"
QUERY="source ~/.bashrc; bcftools query -f '%CHROM %POS %REF %ALT %ID\n'"
COMMAND=paste(QUERY, INPUT)
df0 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt", "id"))
df0[,1:2] %>%
  filter(duplicated(.) | duplicated(., fromLast = TRUE)) %>%
  distinct() -> output
# write to a tab delimited file
OUTPUT=args[2]
fwrite(output, file=OUTPUT, quote=FALSE, sep="\t",
       row.names = FALSE, col.names = FALSE)