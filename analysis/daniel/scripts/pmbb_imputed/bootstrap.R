library(boot)

#Usage:
#Rscript scripts/bootstrap.R <.best from PRSice> <outdir>

#https://www.statmethods.net/advstats/bootstrapping.html


#parse args and format
args = commandArgs(trailingOnly=TRUE)

pref <- args[1]
outdir <- args[2]

covs <- read.table("covs.txt", head=T)
hl <- read.table("deghl.txt", head=T)
scores <- read.table("prscs/scores.txt", head=T)
anc <- read.table("IID_Anc.txt", head=T)
snhl <- read.table("snhl_plus1.txt", head=T)

merged <- merge(covs, hl)
merged <- merge(merged, snhl)
merged <- merge(merged, anc)
merged <- merge(merged, scores)

eur <- subset(merged, Ancestry == "EUR")
afr <- subset(merged, Ancestry == "AFR")

# function to obtain PRS R-Squared from the data
rsq <- function(data, indices) {
	d <- data[indices,] # allows boot to select sample
	
	fit <- lm(DegreeHL ~ Age + AgeSq + Gender, data=d)
	rsq_covs <- summary(fit)$r.square

	fit <- lm(DegreeHL ~ PRS + Age + AgeSq + Gender, data=d)
	rsq_full <- summary(fit)$r.square

	rsq_prs <- rsq_full - rsq_covs

	return(rsq_prs)
}


# bootstrapping with 5000 replications
results <- boot(data=eur, statistic=rsq, R=5000)


# view results
sink(paste(outdir, pref, "_summary.txt", sep=""))
results
cat("mean:")
mean(mean(results$t))
sink()

write.table(results$t, paste(outdir, pref, "_boot.txt", sep=""))

png(filename = paste(outdir, pref, ".png", sep=""))
plot(results)
dev.off()
