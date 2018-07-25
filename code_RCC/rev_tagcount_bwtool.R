#!/usr/bin/env Rscript

## R script for flipping the tagcounts generated from bwtool for motifs on the reverse (minus) strand

args<-commandArgs(trailingOnly=T);
fwd_name = args[1]
rev_name = args[2]
sites_name = args[3]

cat("Flipping tagcounts for motifs on - strand ... \n")

fwd_count.df <- read.table(fwd_name)
rev_count.df <- read.table(rev_name)
sites.df <- read.table(sites_name)

if(sum(sites.df[,2] != fwd_count.df[,2]) != 0){

  stop("mismatch sites!")

} else{

  # Extract the count values
  fwd_count.m <- fwd_count.df[, 6:ncol(fwd_count.df)]
  rev_count.m <- rev_count.df[, 6:ncol(rev_count.df)]

  fwd_output <- fwd_count.m
  rev_output <- rev_count.m

  # For motifs match to the minus strand, flip the fwd and rev tagcounts, and reverse the counts
  idx_minusStrand <- which(sites.df[,6] == "-")

  fwd_output[idx_minusStrand, ] <- t(apply(rev_count.m[idx_minusStrand, ], 1, rev))
  rev_output[idx_minusStrand, ] <- t(apply(fwd_count.m[idx_minusStrand, ], 1, rev))

  # Add the sites info to the first few columns
  fwd_count.df <- cbind(sites.df[,c(1:3,6)], fwd_output)
  rev_count.df <- cbind(sites.df[,c(1:3,6)], rev_output)

  write.table(fwd_count.df, fwd_name, append = F, quote = F, sep = " ", row.names = F, col.names = F)
  write.table(rev_count.df, rev_name, append = F, quote = F, sep = " ", row.names = F, col.names = F)

  cat("Finished flipping tagcounts for motifs on - strand ... \n")

}
