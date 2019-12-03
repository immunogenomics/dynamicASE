wd <- "/Users/mgutierr/Documents/work/rnaseq/cd4timelinepilot/DQB1_LS/rnaseq/expression/"
setwd(wd)
library(preprocessCore)

# load data with expression matrix and meta data matrix
load("/Users/mgutierr/Documents/work/rnaseq/cd4timelinepilot/DQB1_LS/rnaseq/expression/cd4tcLS_rnaseq_expr_swap_filt.rda")

# data can also be read from text matrices from supplementary files in GEO: GSE140244
# we use the "normalized expression" levels to plot here using matrix log2tpm
# these are log2(tpm+1) values

dim(log2tpm)

### filter expression matrix to include onlye "expressed genes"

table(m$timePoint)
thres <- min(table(m$timePoint))
filter <- apply(log2tpm, 1, function(x) length(x[x>2])>=thres) 
filter <- apply(log2tpm, 1, function(x) length(x[x>1])>=thres) 
filtered <- log2tpm[filter,]
genes <- rownames(filtered)[grep("ENS", rownames(filtered))]
length(genes)

log2tpm <- filtered


######### FUNCTION TO PLOT BOXPLOTS

PLOT_GENE <- function(gene){
  
  reorder <- c(which(m$timePoint %in% "0"), which(m$timePoint %in% "2"), which(m$timePoint %in% "4"), which(m$timePoint %in% "8"), which(m$timePoint %in% "12"), which(m$timePoint %in% "24"), which(m$timePoint %in% "48"), which(m$timePoint %in% "72"))
  m <- m[reorder,]
  log2tpm <- log2tpm[,reorder]
  
  cbPalette <- "olivedrab3"
  
  #geneID <- rownames(match)[which(match$GENE_NAME %in% gene)]
  
  x <- as.numeric(log2tpm[geneID,])
  yrange <- c(0, max(x))
  
  boxplot(list(t0=x[which(m$timePoint %in% "0")], t2=x[which(m$timePoint %in% "2")], t4=x[which(m$timePoint %in% "4")], t8=x[which(m$timePoint %in% "8")], t12=x[which(m$timePoint %in% "12")], t24=x[which(m$timePoint %in% "24")], t48=x[which(m$timePoint %in% "48")], t72=x[which(m$timePoint %in% "72")]), 
          col = cbPalette, ylab = "log2(tpm+1)", cex.names = 1.5, cex.lab = 1.5, 
          cex.axis = 1.5, main = gene, cex.main = 1.5, ylim = yrange)
}

####### GENE EXAMPLES TO PLOT

x11(width = 6, height = 5)
PLOT_GENE("IL10")
dev.print(dev = pdf, file = "boxplots_IL10.pdf")

x11(width = 6, height = 5)
PLOT_GENE("UBASH3A")
dev.print(dev = pdf, file = "boxplots_UBASH3A.pdf")

