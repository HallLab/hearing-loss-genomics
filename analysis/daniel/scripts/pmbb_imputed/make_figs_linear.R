library(qqman)
linear <- read.table("/project/ritchie07/personal/daniel/HearingLoss/PMBB_Imputed/results/qq_formatted/allChr_linear_maf.01.txt", head=T, sep="\t")

png("/project/ritchie07/personal/daniel/HearingLoss/PMBB_Imputed/results/qq_formatted/plots/linear_maf.01.png", width = 600, height = 600)
qq(linear$P, main="Degree HL, MAF > .01")
dev.off()
