#!/usr/bin/env Rscript

library(CENTIPEDE)

args <- commandArgs(trailingOnly=T);
tf_name <- args[1]
pwm_id <- args[2]
bam_name <- args[3]
thresh_mapability <- as.numeric(args[4]) # mapability cutoff
thresh_PWMscore <- as.numeric(args[5]) # PWM score cutoff
thresh_pvalue <- args[6] # FIMO motif match p-value threshold

pwm_name <- paste(tf_name, pwm_id, thresh_pvalue, sep = "_")

dir_bam <- "/project/mstephens/ATAC_DNase/ATAC-seq_Olivia_Gray/ATAC-seq_BAMfiles/"
dir_sites <- paste0("/project/mstephens/ATAC_DNase/motif_sites_JASPAR2018/candidate_sites/", thresh_pvalue)
dir_count_matrix <- "/project/mstephens/ATAC_DNase/ATAC-seq_Olivia_Gray/results/ATAC-seq_count_matrix/"
dir_predictions <- "/project/mstephens/ATAC_DNase/ATAC-seq_Olivia_Gray/results/centipede_predictions/"

flank <- 100

##### Functions #####
## load and combine count matrices
load_combine_counts <- function(bam_basename, pwm_name, dir_count_matrix){
  cat("Loading count matrices ... \n")
  counts_fwd.df <- read.table(paste0(dir_count_matrix, "/", pwm_name, "/", pwm_name, "_", bam_basename, "_fwdcounts.m.gz"))
  counts_rev.df <- read.table(paste0(dir_count_matrix, "/", pwm_name, "/", pwm_name, "_", bam_basename, "_revcounts.m.gz"))

  ## the first 5 columns from "bwtool extract" are chr, start, end, name, and the number of data points
  counts_fwd.df <- counts_fwd.df[, -c(1:5)]
  counts_rev.df <- counts_rev.df[, -c(1:5)]

  colnames(counts_fwd.df) <- paste0("fwd", 1:ncol(counts_fwd.df))
  colnames(counts_rev.df) <- paste0("rev", 1:ncol(counts_rev.df))

  counts_combined.m <- as.matrix(cbind(counts_fwd.df, counts_rev.df))

  return(counts_combined.m)
}

## select candidate sites by mapability and PWM score cutoffs
select_sites <- function(sites.df, thresh_mapability=NULL, thresh_PWMscore=NULL, readstats_name=NULL){
  #  cat("loading sites ...\n")
  cat("Select candidate sites with mapability >=", thresh_mapability, ", PWM >=", thresh_PWMscore, "\n")

  if(!is.null(thresh_mapability)){
    idx_mapability <- (sites.df[,"mapability"] >= thresh_mapability)
  }else{
    idx_mapability <- rep(TRUE, nrow(sites.df))
  }

  if(!is.null(thresh_PWMscore)){
    idx_pwm <- (sites.df[,"pwm_score"] >= thresh_PWMscore)
  }else{
    idx_pwm <- rep(TRUE, nrow(sites.df))
  }

  if(!is.null(readstats_name)){
    readstats.df <- read.table(readstats_name, header = F)
    ## if the readstats.df contains chrY, then it means the cell type is male, then the candidate sites should contain chrY,
    ## otherwise, the cell type is female, then the candidate sites on chrY should be removed.
    if( "chrY" %in% readstats.df[,1] ){
      cat("include chrY sites \n")
      idx_chr <- (sites.df[,1] != "")
    }else{
      cat("chrY NOT in the bam file, filter out chrY sites \n")
      ## remove chrY from candidate (motif) sites
      idx_chr <- (sites.df[,1] != "chrY")
    }

  }else{
    idx_chr <- rep(TRUE, nrow(sites.df))
  }

  idx_select <- idx_mapability & idx_pwm & idx_chr

  return(idx_select)
}


## plot data and prediction heatmap
heatmap_data_results <- function(pwm, data, results.df, rank, tf_name, title_name, label_names, zMax_data, zMax_results){

  data[data > zMax_data] <- zMax_data

  results_truncated.df <- results.df
  results_truncated.df[results_truncated.df > zMax_results] <- zMax_results

  myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,0.5)))(1000)

  colNumber = 2 + ncol(results.df)

  layout(matrix(c(1:colNumber), nrow = 1, ncol = colNumber, byrow = F), widths=c(1.4, 4, rep(1.4, ncol(results.df))), heights = rep(1,colNumber))
  cex.text = 0.7
  padj.x = -1.6
  padj.y = 1.2
  par(oma = c(2,1,1,1))

  cat('plotting PWM ... \n')
  par(mai=c(0.1,0.2,0.5,0.2))
  image(t(as.matrix(pwm[rank])), col = myColors, axes=FALSE)
  box()
  mtext(label_names[1], 3, line=0.5, cex =cex.text)
  mtext(paste(length(pwm), tf_name, 'candidate sites'), 2, line= 1, cex = cex.text)

  cat('plotting data ... \n')
  par(mai=c(0.1,0.2,0.5,0.2))
  x = c(1:ncol(data))
  y = c(1:nrow(data))
  image(x, y , z = t(data[rank,]), col = myColors, axes=FALSE, zlim = c(0,zMax_data))
  box()
  axis(1,at=c(1, 101, ncol(data)-100, ncol(data)),labels=c(-100, '','' ,100),padj=padj.x, cex.axis= cex.text, tck=-0.03, tick = T, cex = cex.text)
  mtext(label_names[2], 3,line=0.5,cex =cex.text)
  mtext('Dist. to motif (bp)', 1, line=1.2, cex = cex.text-0.1)

  cat('plotting results ... \n')
  for ( i in 1: ncol(results_truncated.df)){
    par(mai=c(0.1,0.2,0.5,0.2))
    image(t(as.matrix(results_truncated.df[rank,i])), col = myColors, axes=FALSE)
    box()
    mtext(label_names[2+i],3,line=0.5,cex = cex.text)
  }

  mtext(title_name, side=3, outer=TRUE, line= -0.5, cex = cex.text + 0.1)

}

###### begins here ######
bam_basename <- tools::file_path_sans_ext(basename(bam_name))

cat("bam file:", bam_basename, "\n")

cat("pwm_name:", pwm_name, "\n")

filename_sites <- paste0(dir_sites, "/", pwm_name, "_flank", flank, "_fimo_sites.bed")

sites.df <- read.table(filename_sites, header = T, comment.char = "!", stringsAsFactors = F)
colnames(sites.df) = c("chr", "start", "end", "name", "pwm_score", "strand", "p_value", "mapability")

readstats_name <-  paste0(dir_bam, "/", bam_basename, ".idxstats.txt")

idx_select <- select_sites(sites.df, thresh_mapability, thresh_PWMscore, readstats_name)

sites.df <- sites.df[idx_select, c("chr", "start", "end", "name", "pwm_score", "strand", "p_value")]
cat("Number of sites:", nrow(sites.df), "\n")

counts_combined.m <- load_combine_counts(bam_basename, pwm_name, dir_count_matrix)
counts_combined.m <- counts_combined.m[idx_select,]
cat("Dimension of", dim(counts_combined.m), "\n")

if(nrow(counts_combined.m) != nrow(sites.df)){
  stop("Sites not matched!")
}

centFit <- fitCentipede(Xlist = list(Cuts=as.matrix(counts_combined.m)),
                        Y = as.matrix(cbind(IntCept=1, PWM = sites.df$pwm_score)))

site_predictions.df <- data.frame(sites.df,
                                  CentPostPr = signif(centFit$PostPr, 6),
                                  CentLogRatios = signif(centFit$LogRatios, 6))

# logOdds <- log(centFit$PostPr/(1-centFit$PostPr))

dir.create(paste0(paste0(dir_predictions, "/", pwm_name)), recursive = T, showWarnings = F)

write.table(site_predictions.df, paste0(dir_predictions, "/", pwm_name, "/", pwm_name, "_", bam_basename, "_predictions.txt"),
            col.names = T, row.names = F, quote = F)

pdf(paste0(dir_predictions, "/", pwm_name, "/", pwm_name, "_", bam_basename, "_centipede_figures.pdf"))
imageCutSites(counts_combined.m[order(centFit$PostPr),][c(1:100, (dim(counts_combined.m)[1]-100):(dim(counts_combined.m)[1])),])
plotProfile(centFit$LambdaParList[[1]], Mlen = ncol(counts_combined.m)/2 - 200)
hist(centFit$PostPr)
hist(centFit$LogRatios)
dev.off()

## heatmap
condition <- unlist(strsplit(bam_basename, split = "_"))[1]

counts.m <- as.matrix(counts_combined.m[, 1:(ncol(counts_combined.m)/2)] +
                        counts_combined.m[, (ncol(counts_combined.m)/2+1):ncol(counts_combined.m)])

results.df <- data.frame(PostPr = centFit$PostPr, bound = as.integer(centFit$PostPr > 0.99))
label_names <- c('PWM\nscore', paste0('ATAC-seq cuts\n', condition), paste0('P(bound)\n', condition), paste0('Bound label\n', condition))

filename_heatmap <- paste0(dir_predictions, "/", pwm_name, "/", pwm_name, "_", bam_basename, "_centipede_heatmap_sortby_postpr.png")
png(filename_heatmap, width = 4, height = 3, units = 'in', res = 300)

heatmap_data_results(sites.df$pwm_score, counts.m, results.df, order(centFit$PostPr),
                    tf_name, title_name = "", label_names, zMax_data = 5, zMax_results = 1)
dev.off()


cat("Finished CENTIPEDE prediction for: ", pwm_name, "in", bam_name, "\n")
cat("Results saved at", dir_predictions, "\n")

