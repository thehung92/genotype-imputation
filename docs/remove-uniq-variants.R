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
# args <- c("/Users/hung/Data/DGV4VN/Impute_protocol/run_4/data/test920.chr20_QC1.vcf.gz", "/Users/hung/Data/1KGP_LC/panel/1KGP3.chr20_QC1.vcf.gz", "./")
#
library(data.table)
library(tidyverse)
# INPUT1="/Users/hung/Data/DGV4VN/Impute_protocol/run_4/data/test920.chr20_QC1.vcf.gz"
INPUT1=args[1]
QUERY="source ~/.bashrc; bcftools query -i 'INFO/AC=1' -f '%CHROM %POS %REF %ALT\n'"
COMMAND=paste(QUERY, INPUT1)
df1 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt"))
df1 <- df1 %>%
  mutate(id=paste(chrom, pos, ref, alt, sep=":"))
# INPUT2 from 1KGP panel
INPUT2=args[2]
QUERY="source ~/.bashrc; bcftools query -i 'INFO/AC=1' -f '%CHROM %POS %REF %ALT\n'"
COMMAND=paste(QUERY, INPUT2)
df2 <- fread(cmd=COMMAND, sep=" ", header=FALSE,
             col.names=c("chrom", "pos", "ref", "alt"))
df2 <- df2 %>%
  mutate(id=paste(chrom, pos, ref, alt, sep=":"))
# take uniq variant with AC=1 from df1
as_tibble(df1) %>% mutate(RA=paste(ref,alt, sep=",")) %>%
  filter(! id %in% df2$id) %>%
  select(1,2) -> output1
# take uniq variant with AC=1 from df2
as_tibble(df2) %>% mutate(RA=paste(ref,alt, sep=",")) %>%
  filter(! id %in% df1$id) %>%
  select(1,2) -> output2
# write output and compressed and tabix index for bcftools
OUTDIR=args[3]
#
OUTFILE <- basename(INPUT1) %>% gsub("_.*","_uniqAC1.txt", x=.)
OUTPUT1=paste0(OUTDIR, OUTFILE)
fwrite(output1, file=OUTPUT1, quote=FALSE, sep="\t",
       row.names = FALSE, col.names = FALSE)
# COMMAND <- paste("bgzip", OUTPUT1, "; tabix -s1 -b2 -e2", paste0(OUTPUT1,".gz"))
# system(COMMAND)
#
OUTFILE <- basename(INPUT2) %>% gsub("_.*","_uniqAC1.txt", x=.)
OUTPUT2=paste0(OUTDIR, OUTFILE)
fwrite(output2, file=OUTPUT2, quote=FALSE, sep="\t",
       row.names = FALSE, col.names = FALSE)
# COMMAND <- paste("bgzip", OUTPUT2, "; tabix -s1 -b2 -e2", paste0(OUTPUT2,".gz"))
# system(COMMAND)
#
print("Done!")
