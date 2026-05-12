import sys

genes = []
with open(sys.argv[1]) as fp: #allGenes/HL_caseAudAndPhecode/meta_format/chr22.txt
	genes = fp.readline().rstrip().split(",")

toprint = "Population"
for i in range(1, len(genes)):
	toprint += "\t" + genes[i] + "_beta\t" + genes[i] + "_SE\t" + genes[i] + "_p"
print(toprint)

eur_effects = "EUR"
with open(sys.argv[2]) as fp: #e.g., allGenes/HL_rmAudNA/results_to_meta/rmAudNA_EUR_chr12.txt
	for line in fp:
		if line.startswith("eur[,"):
			line = line.rstrip().split()
			beta = line[2]
			se = line[3]
			p = line[5]
			eur_effects += "\t" + beta + "\t" + se + "\t" + p
print(eur_effects)
	
afr_effects = "AFR"
with open(sys.argv[3]) as fp: #e.g., allGenes/HL_rmAudNA/results_to_meta/rmAudNA_AFR_chr12.txt
	for line in fp:
		if line.startswith("afr[,"):
			line = line.rstrip().split()
			beta = line[2]
			se = line[3]
			p = line[5]
			afr_effects += "\t" + beta + "\t" + se + "\t" + p
print(afr_effects)
