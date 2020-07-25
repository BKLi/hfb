## -------------------------------------------------------------

# Created on 12-30-2019
# @ author: bingkun
# @ project HMAMGA

# liftover start/end of interactions seperately
# add ID to interactions
# S1 of pipeline

# previous step: Sep end1/end2 of interaction using awk
## -------------------------------------------------------------


library(rtracklayer)
library(GenomicRanges)
library(R.utils)


cellType <- 'RGC'


infile_all <- sprintf("C:\\Users\\libin\\UCSF\\MAGMA\\%s.5k.downsample.bedpe", cellType)
infile_1 <- sprintf("C:\\Users\\libin\\UCSF\\MAGMA\\%s.5k.downsample.bedpe.end1", cellType)
infile_2 <- sprintf("C:\\Users\\libin\\UCSF\\MAGMA\\%s.5k.downsample.bedpe.end2", cellType)

outfile_all <- sprintf("C:\\Users\\libin\\UCSF\\MAGMA\\%s.5k.downsample.bedpe.withID", cellType)
outfile_1 <- sprintf("C:\\Users\\libin\\UCSF\\MAGMA\\%s.5k.downsample.bedpe.hg19.end1", cellType)
outfile_2 <- sprintf("C:\\Users\\libin\\UCSF\\MAGMA\\%s.5k.downsample.bedpe.hg19.end2", cellType)

##------------------------------------------------------------------
infile_df_all <- read.table(infile_all, sep = "\t", header = TRUE)

infile_df_1 <- read.table(infile_1, sep = "\t", header = FALSE)
names(infile_df_1) <- c("CHR1","START1","END1")

infile_df_2 <- read.table(infile_2, sep = "\t", header = FALSE)
names(infile_df_2) <- c("CHR2","START2","END2")

id_all <- paste0(cellType, rownames(infile_df_all))
id1 <- paste0(cellType, rownames(infile_df_1))
id2 <- paste0(cellType, rownames(infile_df_2))

infile_df_all <- cbind(ID=id_all, infile_df_all)
infile_df_1 <- cbind(ID=id1, infile_df_1)
infile_df_2 <- cbind(ID=id2, infile_df_2)

# convert to Grange object required for liftOver
infile_grange_1 <- makeGRangesFromDataFrame(infile_df_1, keep.extra.columns = TRUE, ignore.strand = TRUE, seqnames.field = "CHR1",start.field = "START1",end.field = "END1" )
infile_grange_2 <- makeGRangesFromDataFrame(infile_df_2, keep.extra.columns = TRUE, ignore.strand = TRUE, seqnames.field = "CHR2",start.field = "START2",end.field = "END2" )

##------------------------------------------------------------------
# liftOver end1
chainpath = getwd()
chainfile <- file.path(chainpath, "hg38ToHg19.over.chain")
chainfile <- import.chain(chainfile)

seqlevelsStyle(infile_grange_1) = "UCSC"
out.coords.1 <- unlist(liftOver(infile_grange_1,chainfile))
genome(out.coords.1) = "hg19"
out.coords.1 <- data.frame(out.coords.1)
# colnames(out.coords.1) <- c("chr1", "start1", "end1",	"width", "strand", "ID", "chr2", "start2", "end2", "count", "expected", "fdr",	"ClusterLabel",	"ClusterSize",	"ClusterType", "ClusterNegLog10P", "ClusterSummit"
# )
colnames(out.coords.1) <- c("chr1", "start1", "end1",	"width", "strand", "ID")
out.coords.1 <- out.coords.1[, c("chr1", "start1", "end1", "ID")]
write.table(out.coords.1, sep="\t", outfile_1, quote = FALSE, row.names = FALSE)

##-----------------------------------------------------------------
# liftOver end2
chainpath = getwd()
chainfile <- file.path(chainpath, "hg38ToHg19.over.chain")
chainfile <- import.chain(chainfile)

seqlevelsStyle(infile_grange_2) = "UCSC"
out.coords.2 <- unlist(liftOver(infile_grange_2,chainfile))
genome(out.coords.2) = "hg19"
out.coords.2 <- data.frame(out.coords.2)
#colnames(out.coords.2) <- c("chr2", "start2", "end2",	"width", "strand", "ID", "chr1", "start1", "end1", "count", "expected", "fdr",	"ClusterLabel",	"ClusterSize",	"ClusterType", "ClusterNegLog10P", "ClusterSummit"
#)
colnames(out.coords.2) <- c("chr2", "start2", "end2",	"width", "strand", "ID")
out.coords.2 <- out.coords.2[, c("chr2", "start2", "end2", "ID")]
write.table(out.coords.2, sep="\t", outfile_2, quote = FALSE, row.names = FALSE)

# write out original interactions with ID
write.table(infile_df_all, sep="\t", outfile_all, quote = FALSE, row.names = FALSE)