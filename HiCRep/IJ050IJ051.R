library(hicrep)

args <- commandArgs(trailingOnly=TRUE)

print ('res=5000,chr=19,h=20,cutoff=5000000')

IJ050<-read.delim2(args[1],header =FALSE,sep='\t')
IJ051<-read.delim2(args[2],header =FALSE,sep='\t')

IJ050 <- depth.adj(IJ050, 1000000, 5000, out = 0)
IJ051 <- depth.adj(IJ051, 1000000, 5000, out = 0)
IJ050IJ051 <- prep(IJ050,IJ051,5000,20,5000000)
SCC.IJ050IJ051 <- get.scc(IJ050IJ051,5000,5000000)
print ('SCC.IJ050IJ051')
SCC.IJ050IJ051[[3]]
SCC.IJ050IJ051[[4]]

