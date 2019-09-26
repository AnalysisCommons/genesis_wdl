# install/load packages
packages <- c("qqman","data.table","RColorBrewer")
lapply(packages, library, character.only = TRUE)

# make a Quantile Quantile plot, with two types of points
make_qq_frq <- function(pvals1, pvals2, main = "QQ plot"){
  pvals1 <- sort(-log10(pvals1[pvals1 > 0]))
  pvals2 <- sort(-log10(pvals2[pvals2 > 0]))
  ymax <- ceiling(max(max(pvals1), max(pvals2)))
  
  plot(x = qexp(ppoints(length(pvals1)))/log(10), y = pvals1, xlab = "Expected", ylab = "Observed", main = main, col = "#E69F00", cex = .8, bg = "#E69F00", pch = 21, ylim = c(0, ymax))
  abline(0, 1, lty = 2)
  points(x = qexp(ppoints(length(pvals2)))/log(10), y = pvals2, col = "#56B4E9", cex = .8, bg = "#56B4E9", pch = 21)
}

# make a Quantile Quantile plot
make_qq <- function(pvals, main = "QQ plot"){
  pvals <- sort(-log10(pvals[pvals > 0]))
  plot(x = qexp(ppoints(length(pvals)))/log(10), y = pvals, xlab = "Expected", ylab = "Observed", main = main, col = "#000000", cex = .8, bg = "#000000", pch = 21, ylim = c(0, ceiling(max(pvals))))
  abline(0, 1, lty = 2)
}

# make the full summary plot with two QQs and one MH
make_summary_plot <- function(data, pval_col = "pval", alt_frq_col = "freq", chr_col = "chr", pos_col = "pos"){  
  make_qq(data[[pval_col]], main = " ")
  legend('topleft',c(paste0('ALL ',lam.new(data[[pval_col]]))),col=c("#000000"), pch=c(21), bty = 'n')
  make_qq_frq(data[data[[alt_frq_col]] >= 0.01,][[pval_col]], data[data[[alt_frq_col]] < 0.01,][[pval_col]], main = " ")
  legend('topleft', c(paste0('MAF >= 1%  ', lam.new(data[data[[alt_frq_col]] >= 0.01,][[pval_col]])), paste0('MAF < 1%  ', lam.new(data[data[[alt_frq_col]] < 0.01,][[pval_col]]))), col = c("#E69F00", "#56B4E9"), pch = c(21,21), pt.bg = c("#E69F00", "#56B4E9"), , bty = 'n')
  manhattan(data, chr = chr_col, bp = pos_col, p = pval_col, main = "All variants", suggestiveline = -log10(5e-5), genomewideline = -log10(5e-8))
}

# make a small summary plot with one qq and one mh
make_small_summary_plot <- function(data, pval_col = "pval", alt_frq_col = "freq", chr_col = "chr", pos_col = "pos"){
  make_qq(data[[pval_col]], main = " ")
  legend('topleft',c(paste0('ALL ',lam.new(data[[pval_col]]))),col=c("#000000"), pch=c(21), bty = 'n')
  manhattan(data, chr = chr_col, bp = pos_col, p = pval_col, main = "All variants", suggestiveline = -log10(5e-5), genomewideline = -log10(5e-8))
}

# Calculate genomic inflation
lam.new <- function(x,p=.5){
  x = x[!is.na(x)]
  x.quantile <- quantile(x,p)
  round((qchisq(1-x.quantile,1)/qchisq(p,1)),2)
}

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
pval.threshold <- as.numeric(input_args[1])
label <- input_args[2]
assoc.files <- unlist(strsplit(input_args[3], ","))

# Load result files and concat
assoc <- do.call(rbind, lapply(assoc.files, fread, data.table = F, stringsAsFactors = F))

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

# Write out all results
fwrite(assoc, paste0(label, ".all.assoc.csv"), sep=",", row.names = F)

# Write out the top results
top.assoc <- assoc[assoc[,pval] < pval.threshold, ]
if (nrow(top.assoc) == 0){
  fwrite(list(), paste0(label, ".top.assoc.csv"), sep=",", row.names = F)
} else {
  fwrite(top.assoc, paste0(label, ".top.assoc.csv"), sep=",", row.names = F)  
}

# Generate summary plots
if (any(assoc$freq < 0.01)) {
  png(filename = paste0(label,".association.plots.png"), width = 8, height = 8, units = "in", res = 400, type = "cairo")
  layout(matrix(c(1,2,3,3),nrow=2,byrow = T))
  make_summary_plot(assoc, pval_col = "P")
  dev.off()
} else {
  png(filename = paste0(label,".association.plots.png"), width = 12, height = 4, units = "in", res = 400, type = "cairo")
  layout(matrix(c(1,2,2),nrow=1,byrow = T))
  make_small_summary_plot(assoc, pval_col = "P")
  dev.off()
}

