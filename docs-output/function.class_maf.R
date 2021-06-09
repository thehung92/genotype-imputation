class_maf <- function(x, y) {
  # x is the database name, y is the vector of threshold
  # read MAF from 1KGP3
  file=paste0("~/Data/1KGP_LC/frequency/MAF_",x,".chr20.txt")
  df <- fread(file=file, sep="\t", header=TRUE)
  df %>%
    cbind(cut(df$MAF, breaks=y)) %>%
    as_tibble() -> df
  colnames(df)[3] <- "MAF_group"
  return(df)
}
