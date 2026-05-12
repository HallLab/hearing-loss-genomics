import sys

gene = ""
beta = ""
se = ""
p = ""
count = 0
print("Gene\tBeta\tSE\tp")
with open(sys.argv[1]) as fp: #e.g. allGenes/HL_rmAudNA/meta_results/chr8_notFormatted.txt
	for line in fp:
		if line.startswith('[1] "'):
			gene = line.rstrip().replace('[1] "', '')
			gene = gene.replace('_beta"', '')
			#print(gene)
			beta = ""
			se = ""
			p = ""
			count = 0
		elif count == 1:
			beta = line.rstrip().split()[1]
		elif count == 2:
			se = line.rstrip().split()[1]
		elif count == 3:
			p = line.rstrip().split()[1]
			print(gene + "\t" + beta + "\t" + se  + "\t" + p)
		count += 1
