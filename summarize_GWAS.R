# install/load packages
library(data.table)

# make a Quantile Quantile plot, with two types of points
make_qq_frq <- function(pvals1, pvals2, main = "QQ plot"){
  pvals1 <- sort(-log10(pvals1[pvals1 > 0]))
  pvals2 <- sort(-log10(pvals2[pvals2 > 0]))
  ymax <- ceiling(max(max(pvals1), max(pvals2)))
  
  plot(x = qexp(ppoints(length(pvals1)))/log(10), y = pvals1, xlab = "Expected", ylab = "Observed", main = main, col = "#E69F00", cex = .8, bg = "#E69F00", pch = 21, ylim = c(0, ymax))
  abline(0, 1, lty = 2)
  points(x = qexp(ppoints(length(pvals2)))/log(10), y = pvals2, col = "#56B4E9", cex = .8, bg = "#56B4E9", pch = 21)
}

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
make_qq <- function(pvals, main = "QQ plot"){
  pvals <- sort(-log10(pvals[pvals > 0]))
  plot(x = qexp(ppoints(length(pvals)))/log(10), y = pvals, xlab = "Expected", ylab = "Observed", main = main, col = "#000000", cex = .8, bg = "#000000", pch = 21, ylim = c(0, ceiling(max(pvals))))
  abline(0, 1, lty = 2)
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
  if (ymax > 7) abline(h = -log(5e-8), col = 'red', lty = 2)
  if (ymax > 4) abline(h = -log(5e-5), col = 'blue', lty = 2)
}


# make the full summary plot with two QQs and one MH
make_summary_plot <- function(data, pval_col = "pval", alt_frq_col = "freq", chr_col = "chr", pos_col = "pos"){  
  make_qq(data[[pval_col]], main = " ")
  legend(x = 'topleft', y = 'topleft', bquote(lambda == .(lam.new(data[[pval_col]]))),col = "#000000", pch = 21, pt.bg = "#000000", bty = 'n', cex = 1.2)
  
  maf.g <- data[data[[alt_frq_col]] >= 0.01,][[pval_col]]
  maf.l <- data[data[[alt_frq_col]] < 0.01,][[pval_col]]
  make_qq_frq(maf.g, maf.l, main = " ")
  legend(x = 'topleft', 
         y = 'topleft', 
         c(as.expression(bquote(lambda[AF<1] == .(lam.new(maf.l)))), as.expression(bquote(lambda[AF>=1] == .(lam.new(maf.g))))),
         col = c("#E69F00", "#56B4E9"), 
         pch = c(21,21), 
         pt.bg = c("#E69F00", "#56B4E9"), 
         cex = 1.2, 
         bty = 'n')
  
  manhattan(data, chr = chr_col, bp = pos_col, p = pval_col, main = "All variants", suggestiveline = -log10(5e-5), genomewideline = -log10(5e-8))
}

make_summary_plot_v2 <- function(data, pval_col = "pval"){  
  make_qq_v2(data, pval_col)
  make_qq_frq_v2(data, pval_col)
  make_manhattan(data, pval_col)
}

# make a small summary plot with one qq and one mh
make_small_summary_plot <- function(data, pval_col = "pval", alt_frq_col = "freq", chr_col = "chr", pos_col = "pos"){
  make_qq(data[[pval_col]], main = " ")
  # legend('topleft',c(paste0('ALL ',lam.new(data[[pval_col]]))),col=c("#000000"), pch=c(21), bty = 'n')
  legend(x = 'topleft', y = 'topleft', bquote(lambda == .(lam.new(data[[pval_col]]))),col=c("#000000"), pch=c(21), bty = 'n', cex = 1.2)
  manhattan(data, chr = chr_col, bp = pos_col, p = pval_col, main = "All variants", suggestiveline = -log10(5e-5), genomewideline = -log10(5e-8))
}

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
assoc.files <- unlist(strsplit(input_args[3], ","))

# Load result files and concat
assoc <- do.call(rbind, lapply(assoc.files, fread, data.table = F, stringsAsFactors = F))

# stop if files are empty
if (nrow(assoc) == 0){
  # fwrite(list(NA), paste0(results.file, ".all_variants.assoc.csv"), sep=",", row.names = F)
  # fwrite(list(NA), paste0(results.file, ".top_variants.assoc.csv"), sep=",", row.names = F)
  write_table(assoc, paste0(results.file, ".all_variants.assoc.csv"))
  write_table(assoc, paste0(results.file, ".top_variants.assoc.csv"))
  png(filename = paste0(results.file,".association.plots.png"), width = 1, height = 1, units = "in", res = 50, type = "cairo")
  dev.off()
} else {
  # get right pval column
  pval <- names(assoc)[grep("pval",names(assoc))][1]
  
  # Make sure the columns are in the right format
  assoc$chr <- sub("^chr", "", as.character(assoc$chr))
  if (any(assoc$chr == "X")) assoc[assoc$chr == "X", "chr"] <- 23
  if (any(assoc$chr == "Y")) assoc[assoc$chr == "Y", "chr"] <- 24
  if (any(assoc$chr == "M")) assoc[assoc$chr == "M", "chr"] <- 25
  assoc$chr <- as.numeric(as.character(assoc$chr))
  assoc$pos <- as.numeric(as.character(assoc$pos))
  assoc$P <- as.numeric(as.character(assoc[,pval]))
  
  # Write out the top results
  top.assoc <- assoc[assoc[,pval] < pval.threshold, ]
  
  # Write out all results
  write_table(assoc, paste0(results.file, ".all_variants.assoc.csv"))
  write_table(top.assoc, paste0(results.file, ".top_variants.assoc.csv"))
  
  # Generate summary plots
  if (any(assoc$freq < 0.01)) {
    png(filename = paste0(results.file,".association.plots.png"), width = 8, height = 8, units = "in", res = 400, type = "cairo")
    layout(matrix(c(1,2,3,3),nrow=2,byrow = T))
    # make_summary_plot(assoc, pval_col = "P")
    make_summary_plot_v2(assoc, pval_col = "P")
    dev.off()
  } else {
    png(filename = paste0(results.file,".association.plots.png"), width = 12, height = 4, units = "in", res = 400, type = "cairo")
    layout(matrix(c(1,2,2),nrow=1,byrow = T))
    # make_small_summary_plot(assoc, pval_col = "P")
    make_small_summary_plot_v2(assoc, pval_col = "P")
    dev.off()
  }
}



