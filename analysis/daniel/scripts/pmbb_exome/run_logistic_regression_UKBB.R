args = commandArgs(trailingOnly=TRUE)

d <- read.table(args[1], head=T)

sink(args[2])
for (i in 11:length(d)){
	print(colnames(d[i]))
	print(i)
	print(summary(glm(HL ~ PC1 + PC2 + PC3 + PC4 + PC5 + Gender + Age + Age_sq + d[,i], d, family = "binomial")))
}
sink()
