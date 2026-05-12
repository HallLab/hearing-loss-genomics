import sys


gene_pos = {}
with open (sys.argv[1]) as fp: #gene_list_regions.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		gene_pos[line[0]] = line[1] + " " + line[2] + " " + line[3]


gene_counts = {}
used = set() #this indv used for this gene
gene = ""
with open(sys.argv[2]) as fp: #e.g. controls_category1.vcf
	for line in fp:
		if not line.startswith("#"):
			line = line.rstrip().split()

			thischr = line[0]
			thispos = float(line[1])

			for g in gene_pos:
				chrom, start, stop = gene_pos[g].split()
				#new gene
				if g != gene and thischr == chrom and float(start) <= thispos <= float(stop):
					gene = g
					used = set() #reset
			
			for i in range(9, len(line)):
				if "1" in line[i] and i not in used:
					if gene not in gene_counts:
						gene_counts[gene] = 1
					else:
						gene_counts[gene] += 1
					used.add(i)




print("Gene\tN_carriers")
for g in sorted(gene_counts):
	print(g + "\t" + str(gene_counts[g]))
