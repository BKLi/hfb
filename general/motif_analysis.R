library("preprocessCore")
library(reshape2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(gplots)

my_palette <- colorRampPalette(c("#CD0000", "#E68080", "white", "#6495ED"))(n = 199)
color.bar <- function(lut, min, max=-min, nticks, ticks=seq(min, max, len=nticks), labels, title='') {
  
  scale = (length(lut)-1)/(max-min)
  
  # dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=title, main=NULL, cex.main=0.85)
  axis(2, ticks, labels=labels, las=2)
  for (i in 1:(length(lut)-1)) {
    
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    
  }
}

jaspar_breaks <- c(seq(2.6e-79, 2.50e-75, length=50),
                    seq(2.501e-75, 2.50e-25,length=50),
                    seq(2.501e-25, 2.50e-10, length=50),
                    seq(2.501e-10, 2.6e-1, length=50))

jaspar_breaks_raw <- c(seq(1e-304, 1e-300, length=50),
                   seq(1.01e-300, 1e-93,length=50),
                   seq(1.01e-93, 1e-17, length=50),
                   seq(1.01e-17, 1, length=50))

jaspar_merged <- read.table('C:\\Users\\libin\\R_projects\\motif_analysis\\jaspar_merged', sep = "\t", header = TRUE)
rownames(jaspar_merged) = jaspar_merged[,1]
jaspar_merged_mat <- data.matrix(jaspar_merged[,-1])
jaspar_merged_norm <- normalize.quantiles(jaspar_merged_mat,copy=TRUE)
rownames(jaspar_merged_norm) <- jaspar_merged[,1]
colnames(jaspar_merged_norm) <- c("RG", "IPC", "eN", "iN")
# melt_jaspar_norm <- melt(jaspar_merged_norm)
# colnames(melt_jaspar_norm) <- c("Motif", "CellType", "P-value")

pdf('C:\\Users\\libin\\R_projects\\motif_analysis\\jaspar_normalized.pdf', height=8, width=4)
heatmap.2(jaspar_merged_norm, 
          Rowv=NA, Colv=NA, srtCol = 0,
          col = my_palette, breaks=jaspar_breaks, 
          lwid = c(0.1,4), lhei = c(0.1,3.5),
          cexRow=0.7, cexCol=0.9, margins=c(5, 7), 
          cellnote = round(-log10(jaspar_merged_norm), digits = 2), notecol="black",
          dendrogram="none", trace="none",  
          density.info="none",key=FALSE)
dev.off()

pdf('C:\\Users\\libin\\R_projects\\motif_analysis\\jaspar_normalized_colorbar.pdf', height=6, width=3)
color.bar(my_palette, min=0, max=10, nticks=5, 
          labels=c(2.6e-79, 2.5e-75, 2.5e-25, 2.5e-10, 2.6e-1), 
          title="")
dev.off()

pdf('C:\\Users\\libin\\R_projects\\motif_analysis\\jaspar_raw.pdf', height=8, width=4)
heatmap.2(jaspar_merged_mat, 
          Rowv=NA, Colv=NA, srtCol = 0,
          col = my_palette, breaks=jaspar_breaks_raw, 
          lwid = c(0.1,4), lhei = c(0.1,3.5),
          cexRow=0.7, cexCol=0.9, margins=c(5, 7), 
          cellnote = round(-log10(jaspar_merged_mat), digits = 2), notecol="black",
          dendrogram="none", trace="none",  
          density.info="none",key=FALSE)
dev.off()

pdf('C:\\Users\\libin\\R_projects\\motif_analysis\\jaspar_raw_colorbar.pdf', height=6, width=3)
color.bar(my_palette, min=0, max=10, nticks=5, 
          labels=c(1e-304, 1e-300, 1e-93, 1e-17, 1), 
          title="")
dev.off()


homer_breaks <- c(seq(2.6e-66, 2.50e-57,length=50),
                   seq(2.501e-57, 2.50e-21,length=50),
                   seq(2.501e-21, 2.10e-10, length=50),
                   seq(2.101e-10, 7.6e-1, length=50))

homer_breaks_raw <- c(seq(0, 1e-192,length=50),
                  seq(1.01e-192, 1e-77,length=50),
                  seq(1.01e-77, 1e-12, length=50),
                  seq(1.01e-12, 1, length=50))

homer_merged <- read.table('C:\\Users\\libin\\R_projects\\motif_analysis\\homer_merged', sep = "\t", header = TRUE)
rownames(homer_merged) = homer_merged[,1]
homer_merged_mat <- data.matrix(homer_merged[,-1])
homer_merged_norm <- normalize.quantiles(homer_merged_mat,copy=TRUE)
rownames(homer_merged_norm) <- homer_merged[,1]
colnames(homer_merged_norm) <- c("RG", "IPC", "eN", "iN")

pdf('C:\\Users\\libin\\R_projects\\motif_analysis\\homer_normalized.pdf', height=8, width=4)
heatmap.2(homer_merged_norm, 
          Rowv=NA, Colv=NA, srtCol = 0,
          col = my_palette, breaks=homer_breaks, 
          lwid = c(0.1,4), lhei = c(0.1,3.5),
          cexRow=0.7, cexCol=0.9, margins=c(5, 7), 
          dendrogram="none", trace="none",  
          density.info="none",key=FALSE, 
          cellnote = round(-log10(homer_merged_norm), digits = 2), notecol="black")
dev.off()

pdf('C:\\Users\\libin\\R_projects\\motif_analysis\\homer_normalized_colorbar.pdf', height=6, width=3)
color.bar(my_palette, min=0, max=10, nticks=5, 
          labels=c(2.6e-66, 2.5e-57, 2.5e-21, 2.1e-10, 7.6e-1), 
          title="")
dev.off()

pdf('C:\\Users\\libin\\R_projects\\motif_analysis\\homer_raw.pdf', height=8, width=4)
heatmap.2(homer_merged_mat, 
          Rowv=NA, Colv=NA, srtCol = 0,
          col = my_palette, breaks=homer_breaks_raw, 
          lwid = c(0.1,4), lhei = c(0.1,3.5),
          cexRow=0.7, cexCol=0.9, margins=c(5, 7), 
          dendrogram="none", trace="none",  
          density.info="none",key=FALSE, 
          cellnote = round(-log10(homer_merged_mat), digits = 2), notecol="black")
dev.off()

pdf('C:\\Users\\libin\\R_projects\\motif_analysis\\homer_raw_colorbar.pdf', height=6, width=3)
color.bar(my_palette, min=0, max=10, nticks=5, 
          labels=c(0, 1e-192, 1e-77, 1e-12, 1), 
          title="")
dev.off()
