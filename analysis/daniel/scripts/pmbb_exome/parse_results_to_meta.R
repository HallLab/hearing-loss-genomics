args = commandArgs(trailingOnly=TRUE)

chrom <- args[1]
pheno <- args[2]

d <- read.table(paste("allGenes/HL_", pheno, "/meta_format/chr", chrom, ".txt", sep=""), sep=",", head=T)
covs <- read.table(paste("allGenes/HL_", pheno, "/covs_HL_nGenesVar_", pheno, ".txt", sep=""), head=T)
anc <- read.table("anc.txt", head=T)

allcovs <- merge(covs, anc)
merged <- merge(allcovs, d, all = T)
eur <- subset(merged, Ancestry == "EUR")

sink(paste("allGenes/HL_", pheno, "/results_to_meta/", pheno, "_EUR_chr", chrom, ".txt", sep=""))
for (i in 15:length(eur)){
  print(summary(glm(SNHL ~ Sex + Age + AgeSq + PC1 + PC2 + PC3 + PC4 + eur[,i], data = eur, family = "binomial")))
}
sink()

afr <- subset(merged, Ancestry == "AFR")
sink(paste("allGenes/HL_", pheno, "/results_to_meta/", pheno, "_AFR_chr", chrom, ".txt", sep=""))
for (i in 15:length(afr)){
  print(summary(glm(SNHL ~ Sex + Age + AgeSq + PC1 + PC2 + PC3 + PC4 + afr[,i], data = afr, family = "binomial")))
}
sink()

