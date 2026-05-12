args = commandArgs(trailingOnly=TRUE)
library(meta)

d <- read.table(args[1], head=T)
outfile <- args[2]

sink(outfile)
for(i in seq(from=2, to=length(d), by=3)){
	model <- metagen(TE=d[,i], seTE=d[,i+1], sm="OR")
	
	print(colnames(d[i]))
	print(model$TE.fixed)
	print(model$seTE.fixed)
	print(model$pval.fixed)
}
sink()
