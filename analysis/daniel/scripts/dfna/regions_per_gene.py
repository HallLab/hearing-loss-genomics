import sys

gene_pos = {}
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()

		
		gene = line[6]
		chrom = line[0]
		start = int(line[1])
		end = int(line[2])

		if gene not in gene_pos:
			gene_pos[gene] = {}
			gene_pos[gene]["chrom"] = chrom
			gene_pos[gene]["start"] = start
			gene_pos[gene]["end"] = end
		else:
			if chrom == gene_pos[gene]["chrom"]:
				if start < gene_pos[gene]["start"] :
					gene_pos[gene]["start"] = start
				if end > gene_pos[gene]["end"]:
					gene_pos[gene]["end"] = end

print("Gene\tChromosome\tStart\tEnd")
for g in sorted(gene_pos):
	print(g + "\t" + str(gene_pos[g]["chrom"]) + "\t" + str(gene_pos[g]["start"]) + "\t" + str(gene_pos[g]["end"]))
