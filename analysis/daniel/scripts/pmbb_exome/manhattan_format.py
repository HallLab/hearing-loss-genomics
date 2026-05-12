import sys

gene_chr_pos = {}
with open(sys.argv[1]) as fp: #gene_list_regions.txt
	for line in fp:
		line = line.rstrip().split()
		gene_chr_pos[line[1]] = line[0] + "\t" + line[2]

print("SNP\tCHR\tPOS\tP\tDirection") #e.g., allGenes/HL_rmAudNA/meta_results_degHL/all_chrom_meta.txt
with open(sys.argv[2]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] in gene_chr_pos:
			if float(line[1]) < 0:
				print(line[0] + "\t" + gene_chr_pos[line[0]] + "\t" + line[3] + "\tProtective")
			elif float(line[1]) >= 0:
				print(line[0] + "\t" + gene_chr_pos[line[0]] + "\t" + line[3] + "\tDeleterious")
			else:
				print(line[0] + "\t" + gene_chr_pos[line[0]] + "\t" + line[3] + "\tNA")
