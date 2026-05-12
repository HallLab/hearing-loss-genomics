args = commandArgs(trailingOnly=TRUE)

d <- read.table(args[1], head=T)

sink(args[2])
for (i in 27:length(d)){
	print(colnames(d[i]))
	print(i)
	print(summary(lm(DegHL_rand ~ PC1 + PC2 + PC3 + PC4 + PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20 + Sex + Age + AgeSq + d[,i], d)))
}
sink()
