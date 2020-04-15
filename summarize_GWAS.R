# install/load packages
library(data.table)
library(dplyr)

# fix a unicode error
Sys.setlocale(category = "LC_ALL","C.UTF-8")

# normalize association DF
normalizeAssoc <- function(df, pval, remove = T){
  # Make sure the columns are in the right format
  df$chr <- sub("chr", "", df$chr)
  chr.levels <- c(1:22, "X")
  cur.chrs <- unique(df$chr)
  df$chr <- factor(df$chr, levels = chr.levels[chr.levels %in% cur.chrs])
  df$pos <- as.numeric(as.character(df$pos))
  df$P <- as.numeric(as.character(df[,pval]))
  
  # remove NA, 0, and inf pvalues
  if (remove == T){
    df <- df[!is.na(df$P),]
    df <- df[df$P > 0,]
    df <- df[!is.infinite(df$P),]  
  }
  return(df)
}

normalizeSW <- function(df, pval, remove = T){
  # Make sure the columns are in the right format
  df$chr <- sub("chr", "", df$chr)
  chr.levels <- c(1:22, "X")
  cur.chrs <- unique(df$chr)
  df$chr <- factor(df$chr, levels = chr.levels[chr.levels %in% cur.chrs])
  df$start <- as.numeric(as.character(df$start))
  df$end <- as.numeric(as.character(df$end))
  df$pos <- ceiling(unlist(apply(df %>% select(start, end), 1, median, na.rm = T)))
  df$P <- as.numeric(as.character(df[,pval]))
  
  # remove NA, 0, and inf pvalues
  if (remove == T){
    df <- df[!is.na(df$P),]
    df <- df[df$P > 0,]
    df <- df[!is.infinite(df$P),]  
  }
  return(df)
}

# make a Quantile Quantile plot, with two types of points
make_qq_frq_v2 <- function(data, pval_col){
  qq1 <- ramwas::qqPlotPrepare(sort(data[[pval_col]][data$freq >= 0.01 & data[[pval_col]] > 0]))
  qq1$lambda <-  NULL
  qq1$col <- "#56B4E9"
  qq2 <- ramwas::qqPlotPrepare(sort(data[[pval_col]][data$freq < 0.01 & data[[pval_col]] > 0]))
  qq2$lambda <-  NULL
  qq2$col <- "#E69F00"
  
  if (qq1$ypvs[1] > qq2$ypvs[2]){
    ymax <- ceiling(max(qq1$ypvs[1], qq1$xpvs[1]))
    axistep <- dplyr::case_when(
      ymax > 60 ~ 20,
      ymax > 40 ~ 10,
      ymax > 20 ~ 5,
      ymax > 10 ~ 4,
      TRUE ~ 2
    )
    
    ramwas::qqPlotFast(qq1, makelegend = F, col = qq1$col, cex = 0.5, lwd = 1.0, ci.level = NULL, axistep = axistep, xlab = "Expected", ylab = "Observed", yaxmax = ymax+axistep)
    if (axistep > floor(qq1$xpvs[1])) axis(1, seq(0, qq1$xpvs[1], 2), lwd = 1.0);
    if (axistep > floor(qq1$ypvs[1])) axis(2, seq(0, qq1$ypvs[1], 2), lwd = 1.0);
    ramwas::qqPlotFast(qq2, newplot = F, makelegend = F, col = qq2$col, cex = 0.5, lwd = 1.0, ci.level = NULL)
  } else {
    ymax <- ceiling(max(qq2$ypvs[1], qq2$xpvs[1]))
    axistep <- dplyr::case_when(
      ymax > 60 ~ 20,
      ymax > 40 ~ 10,
      ymax > 20 ~ 5,
      ymax > 10 ~ 4,
      TRUE ~ 2
    )
    ramwas::qqPlotFast(qq2, makelegend = F, col = qq2$col, cex = 0.5, lwd = 1.0, ci.level = NULL, axistep = axistep, xlab = "Expected", ylab = "Observed", yaxmax = ymax+axistep)
    if (axistep > floor(qq2$xpvs[1])) axis(1, seq(0, qq2$xpvs[1], 2), lwd = 1.0);
    if (axistep > floor(qq2$ypvs[1])) axis(2, seq(0, qq2$ypvs[1], 2), lwd = 1.0);
    ramwas::qqPlotFast(qq1, newplot = F, makelegend = F, col = qq1$col, cex = 0.5, lwd = 1.0, ci.level = NULL)
  }
  legend(x = 'topleft', 
         y = 'topleft', 
         c(as.expression(bquote(lambda[AF>=1 *'%'] == .(lam.new(data[[pval_col]][data$freq >= 0.01 & data[[pval_col]] > 0])))), as.expression(bquote(lambda[AF<1 *'%'] == .(lam.new(data[[pval_col]][data$freq < 0.01 & data[[pval_col]] > 0]))))),
         col = c(qq1$col, qq2$col), 
         pch = c(21,21), 
         pt.bg = c(qq1$col, qq2$col),
         cex = 1.2, 
         bty = 'n')
}

# make a Quantile Quantile plot
make_qq_v2 <- function(data, pval_col){
  qq <- ramwas::qqPlotPrepare(sort(data[[pval_col]][data[[pval_col]] > 0]))
  qq$lambda <-  NULL
  qq$col <- "#000000"
  ymax <- ceiling(max(qq$ypvs[1], qq$xpvs[1]))
  axistep <- dplyr::case_when(
    ymax > 60 ~ 20,
    ymax > 40 ~ 10,
    ymax > 20 ~ 5,
    ymax > 10 ~ 4,
    TRUE ~ 2
  )
  ramwas::qqPlotFast(qq, makelegend = F, col = qq$col, cex = 0.5, lwd = 1.0, ci.level = NULL, axistep = axistep, xlab = "Expected", ylab = "Observed", yaxmax = ymax+axistep)
  if (axistep > floor(qq$xpvs[1])) axis(1, seq(0, qq$xpvs[1], 2), lwd = 1.0);
  if (axistep > floor(qq$ypvs[1])) axis(2, seq(0, qq$ypvs[1], 2), lwd = 1.0);
  legend(x = 'topleft',
         y = 'topleft', 
         bquote(lambda == .(lam.new(data[[pval_col]]))),
         col = "#000000", 
         pch = 21, 
         pt.bg = "#000000", 
         bty = 'n', 
         cex = 1.2)
}

## make fast mh plot
make_manhattan <- function(data, pval_col){
  man <- ramwas::manPlotPrepare(
    data[[pval_col]],
    data$chr,
    data$pos)
  ymax <- ceiling(max(man$y))
  axistep <- dplyr::case_when(
    ymax > 60 ~ 20,
    ymax > 40 ~ 10,
    ymax > 20 ~ 5,
    ymax > 10 ~ 4,
    TRUE ~ 2
  )
  ramwas::manPlotFast(man, colorSet = c('#474747',"#8f8f8f"), lwd = 1, axistep = axistep, cex = 0.8, yaxmax = ymax+axistep)
  if (ymax > 7) abline(h = -log10(5e-8), col = 'red', lty = 2)
  if (ymax > 4) abline(h = -log10(5e-5), col = 'blue', lty = 2)
}

# make the full summary plot with two QQs and one MH
make_summary_plot_v2 <- function(data, pval_col = "pval"){  
  make_qq_v2(data, pval_col)
  make_qq_frq_v2(data, pval_col)
  make_manhattan(data, pval_col)
}

# make a small summary plot with one qq and one mh
make_small_summary_plot_v2 <- function(data, pval_col = "pval"){  
  make_qq_v2(data, pval_col)
  make_manhattan(data, pval_col)
}

# Calculate genomic inflation
lam.new <- function(x,p=.5){
  x = x[!is.na(x)]
  x.quantile <- quantile(x,p)
  round((qchisq(1-x.quantile,1)/qchisq(p,1)),2)
}

# Helper for writing tables
write_table <- function(data, out.file){
  out.gz <- paste0(out.file, ".gz")
  if (nrow(data) == 0){
    fwrite(list(NA), out.file, sep=",", row.names = F)
    R.utils::gzip(out.file, destname = out.gz, overwrite = T)
  } else {
    fwrite(data, out.file, sep=",", row.names = F)
    R.utils::gzip(out.file, destname = out.gz, overwrite = T)
  }
}

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
pval.threshold <- as.numeric(input_args[1])
results.file <- input_args[2]
test.type <- input_args[3]
agg.file <- input_args[4]
assoc.files <- unlist(strsplit(input_args[5], ","))

# output files
assoc.out_file <- paste0(results.file, ".all_variants.assoc.csv")
top_assoc.out_file <- paste0(results.file, ".top_variants.assoc.csv")
plot.out_file <- paste0(results.file,".association.plots.png")

# Load result files and concat
assoc <- do.call(plyr::rbind.fill, lapply(assoc.files, fread, data.table = F, stringsAsFactors = F))

# stop if files are empty
if (nrow(assoc) == 0){
  write_table(assoc, assoc.out_file)
  write_table(assoc, top_assoc.out_file)
  png(filename = plot.out_file, width = 1, height = 1, units = "in", res = 50, type = "cairo")
  dev.off()
} else {
  # get right pval column
  pval <- names(assoc)[grep("pval",names(assoc))][1]
  
  # single variant tests:
  if (test.type == "Single"){
    # Make sure the columns are in the right format
    assoc.norm <- normalizeAssoc(assoc, pval, remove = T)

    if (nrow(assoc.norm) == 0){
      write_table(assoc, assoc.out_file)
      write_table(assoc, top_assoc.out_file)
      png(filename = plot.out_file, width = 1, height = 1, units = "in", res = 50, type = "cairo")
      dev.off()
    } else {
      # Write out the top results
      top.assoc <- assoc[assoc[,pval] < pval.threshold, ]
      
      # Write out all results
      write_table(assoc, assoc.out_file)
      write_table(top.assoc, top_assoc.out_file)
      
      # Generate summary plots
      if (any(assoc.norm$freq < 0.01)) {
        png(filename = plot.out_file, width = 8, height = 8, units = "in", res = 400, type = "cairo")
        layout(matrix(c(1,2,3,3),nrow=2,byrow = T))
        make_summary_plot_v2(assoc.norm, pval_col = "P")
        dev.off()
      } else {
        png(filename = plot.out_file, width = 12, height = 4, units = "in", res = 400, type = "cairo")
        layout(matrix(c(1,2,2),nrow=1,byrow = T))
        make_small_summary_plot_v2(assoc.norm, pval_col = "P")
        dev.off()
      }
    }
  } else if (agg.file != "NA"){ # tests using an aggregation file
    # load aggregation units
    agg <- fread(agg.file, data.table = F, stringsAsFactors = F) %>%
      rename_all(tolower) %>%
      group_by(group_id) %>%
      mutate(min_pos = min(pos, na.rm = T),
        max_pos = max(pos, na.rm = T)) %>%
        mutate(pos = median(c(min_pos, max_pos)))
    names(agg)[names(agg) %in%  c('chr','chrom','chromosome')] = 'chr' 
    
    # get aggregation format correct
    agg <- agg %>%
      select(group_id, chr, pos, min_pos, max_pos) %>%
      distinct()
  
    # merge with association results
    assoc.norm <- merge(assoc, agg, by.x = 'V1', by.y = 'group_id', all.x = T) %>%
      normalizeAssoc(pval = pval, remove = T)
    names(assoc.norm)[names(assoc.norm) == "V1"] <- "group_id"
    
    if (nrow(assoc.norm) == 0){
      write_table(assoc, assoc.out_file)
      write_table(assoc, top_assoc.out_file)
      png(filename = plot.out_file, width = 1, height = 1, units = "in", res = 50, type = "cairo")
      dev.off()
    } else {
      # write out results
      top.assoc <- assoc[assoc[,pval] < pval.threshold, ]
      write_table(assoc, assoc.out_file)
      write_table(top.assoc, top_assoc.out_file)  
      
      # make plots
      png(filename = plot.out_file, width = 12, height = 4, units = "in", res = 400, type = "cairo")
      layout(matrix(c(1,2,2),nrow=1,byrow = T))
      make_small_summary_plot_v2(assoc.norm, pval_col = "P")
      dev.off()
    }
  } else if (test.type %in% c("Burden", "SKAT", "fastSKAT", "SMMAT",  "SKATO")) { # sliding window
    # normalize assoc
    assoc.norm <- normalizeSW(assoc, pval, remove = T)
    
    if (nrow(assoc.norm) == 0){
      write_table(assoc, assoc.out_file)
      write_table(assoc, top_assoc.out_file)
      png(filename = plot.out_file, width = 1, height = 1, units = "in", res = 50, type = "cairo")
      dev.off()
    } else {
      assoc.norm <- assoc.norm %>%
        arrange(factor(chr, levels = c(as.character(1:22), "X", "Y", "M")), as.numeric(pos))
      
      # get top associations
      top.assoc <- assoc[assoc[,pval] < pval.threshold, ]
      
      # write out results
      write_table(assoc, assoc.out_file)
      write_table(top.assoc, top_assoc.out_file)  
      
      # make plots
      png(filename = plot.out_file, width = 12, height = 4, units = "in", res = 400, type = "cairo")
      layout(matrix(c(1,2,2),nrow=1,byrow = T))
      make_small_summary_plot_v2(assoc.norm, pval_col = "P")
      dev.off()
    }
  } else {
    write_table(assoc, assoc.out_file)
    write_table(assoc, top_assoc.out_file)
    png(filename = plot.out_file, width = 1, height = 1, units = "in", res = 50, type = "cairo")
    dev.off()
  }
}



