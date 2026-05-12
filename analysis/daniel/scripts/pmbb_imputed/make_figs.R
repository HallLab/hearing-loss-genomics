library(qqman)
linear_rand <- read.table("/project/ritchie07/personal/daniel/HearingLoss/PMBB_Imputed/results/qq_formatted/rand_allChr_linear.txt", head=T, sep="\t")

png("/project/ritchie07/personal/daniel/HearingLoss/PMBB_Imputed/results/qq_formatted/plots/linear_rand.png", width = 600, height = 600)
qq(linear_rand$P, main="Linear_rand")
dev.off()

maf.005 <- subset(linear_rand, MAF > .005)
png("/project/ritchie07/personal/daniel/HearingLoss/PMBB_Imputed/results/qq_formatted/plots/linear_rand_maf.005.png", width = 600, height = 600)
qq(maf.005$P, main="Linear_rand, maf.005")
dev.off()

maf.01 <- subset(maf.005, MAF > .01)
png("/project/ritchie07/personal/daniel/HearingLoss/PMBB_Imputed/results/qq_formatted/plots/linear_rand_maf.01.png", width = 600, height = 600)
qq(maf.01$P, main="Linear_rand, maf.01")
dev.off()
