import sys

gene = {}
with open(sys.argv[1]) as fp: #gene_list_regions.txt -- manually checked, none overlap
	for line in fp:
		line = line.rstrip().split()
		gene[line[1]] = {}
		gene[line[1]]["chrom"] = line[0]
		gene[line[1]]["start"] = float(line[2])
		gene[line[1]]["end"] = float(line[3])


id_ind = {}
gene_indv_counts = {}
with open(sys.argv[2]) as fp: #allIndvs_maf.001_noRels_merged.vcf
	for line in fp:
		#get ind of ids
		if line.startswith("#CHROM"):
			line = line.rstrip().split()
			for i in range(9, len(line)):
				id_id[line[0]] = i

	
		#count	
		elif not line.startswith("#"):
			line = line.rstrip().split()

			chrom = line[0]
			pos = float(line[1])

			currgene = "NA"
			for g in gene:
				if chrom == gene[g]["chrom"] and pos >= gene[g]["start"] and pos <= gene[g]["end"]:
					currgene = g

			if currgene not in gene_indv_counts:
				gene_indv_counts = {}
				for indv in id_ind:
					gene_indv_counts[currgene]["indv"] = 0
		
			for i in range(9, len(line)):
				

