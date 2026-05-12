import sys

genes = {}
with open(sys.argv[1]) as fp: #annot_genes_full_funcToInclude.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()

		gene = line[7]
		chrom = line[1]
		start = line[2]
		stop = line[3]	

		if gene not in genes:
			genes[gene] = []
			genes[gene].append(chrom)
			genes[gene].append(start)
			genes[gene].append(stop)
		else:
			if float(start) < float(genes[gene][1]):
				genes[gene][1] = start
			if float(stop) > float(genes[gene][2]):
				genes[gene][2] = stop

for g in sorted(genes):
	print(genes[g][0] + "\t" + g + "\t" + genes[g][1] + "\t" + genes[g][2])
