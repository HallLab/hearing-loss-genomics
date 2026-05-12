library(qqman)
binary <- read.table("/project/ritchie07/personal/daniel/HearingLoss/PMBB_Imputed/results/qq_formatted/allChr_binary_maf.01.txt", head=T, sep="\t")

png("/project/ritchie07/personal/daniel/HearingLoss/PMBB_Imputed/results/qq_formatted/plots/binary_maf.01.png", width = 600, height = 600)
qq(binary$P, main="Binary, maf.01")
dev.off()
