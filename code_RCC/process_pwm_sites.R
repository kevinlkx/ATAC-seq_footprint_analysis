#!/usr/bin/env Rscript

## R script for flanking the motif sites and compute mapability
## use the bigWigAverageOverBed tool from UCSC to compute mapablity
## https://genome.ucsc.edu/goldenpath/help/bigWig.html

args <- commandArgs(trailingOnly=T);
tf_name = args[1]
pwm_id = args[2]
filename_fimo = args[3]
flank = as.integer(args[4])
thresh_pValue = args[5]
dir_sites = args[6]
bw_mapability = args[7] # wgEncodeDukeMapabilityUniqueness20bp.bigWig

cat(dir_sites, '\n')
dir.create(dir_sites, showWarnings = F, recursive = T)

options(scipen=999) ## suppress scientific notations

if(!file.exists(bw_mapability)){
  cat("Warning: Mappability uniqueness file does not exist! \n")
}

cat("processing sites:", filename_fimo, flank, thresh_pValue, "\n")
fimo.df <- read.table(filename_fimo, header = T, stringsAsFactors = F, comment.char = "!", sep = "\t")

if(nrow(fimo.df) < 1){
  stop('No fimo sites matched for', tf_name, pwm_id, '\n')

}

## order sites
chr_order <- paste0("chr", c(1:22, 'X','Y','M'))
fimo.df <- fimo.df[fimo.df$sequence_name %in% chr_order,]
fimo.df$sequence_name <- factor(fimo.df$sequence_name, chr_order, ordered=TRUE)
fimo.df <- fimo.df[order(fimo.df$sequence_name, fimo.df$start, fimo.df$stop),]

## collect fimo output's coordinates, pwm score, strand, p-value, q-value
fimo.bed <- data.frame(chr = fimo.df$sequence_name, start = fimo.df$start, stop = fimo.df$stop,
                       name = paste0("site", c(1:nrow(fimo.df))),
                       score = fimo.df$score, strand = fimo.df$strand,
                       p.value = fimo.df$p.value)

## the fimo output sites are 1-based, so convert it to 0-based coordinates
## make the coordinates the same format as BED, [Start: end), not including the end position, so add 1 to the stop position
fimo.bed$start <- fimo.bed$start -1 - flank
fimo.bed$stop <- fimo.bed$stop -1 + flank + 1

# filter out the sites whose window start before 0 after the 100bp flanking
fimo.bed <- fimo.bed[fimo.bed$start >= 0, ]

## thresholding the sites with p-value
cat("Number of sites after thresh_pValue =", thresh_pValue, ":",nrow(fimo.bed), "-> ")
fimo.bed <- fimo.bed[which(as.numeric(fimo.bed$p.value) <= as.numeric(thresh_pValue)), ]
cat(nrow(fimo.bed), "\n")

filename_sites <- paste0(dir_sites, "/", tf_name, "_", pwm_id, "_", thresh_pValue, "_flank", flank, "_fimo_sites.bed")

if(file.exists(bw_mapability)){
  cat("Compute mapability ... \n")
  ## creating a bed file for the use of mapability
  fimo_mapability.bed <- fimo.bed[,1:6]
  fimo_mapability.bed$name <- paste("site", c(1:nrow(fimo.bed)), sep = "")
  fimo_mapability.bed$score <- 0

  filename_bed_tmp <- paste0(filename_sites,"_bed.tmp")
  write.table(fimo_mapability.bed, filename_bed_tmp, sep = "\t", quote = F, row.names = F, col.names = F)

  filename_mapability_tmp <- paste(filename_sites,"_mapability.tmp", sep = "")
  try(system(paste("bigWigAverageOverBed", bw_mapability, filename_bed_tmp, filename_mapability_tmp)))

  # mean0 column from bigWigAverageOverBed output: average over bases with non-covered bases counting as zeroes
  fimo.bed <- data.frame(fimo.bed, mapability = read.table(filename_mapability_tmp)[,5])

  try(system(paste("rm", filename_bed_tmp)))
  try(system(paste("rm", filename_mapability_tmp)))
}

names(fimo.bed)[1] <- paste0('#', names(fimo.bed)[1])  ## add # to the header
write.table(fimo.bed, file = filename_sites, sep = "\t", quote = F, row.names = F, col.names = T)

cat(paste(nrow(fimo.bed),"sites written in", filename_sites, "\n"))

