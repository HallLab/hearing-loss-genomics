import sys

ind_id = {} #index mapping to ID
id_N_sites = {} #ID mapping to N sites
with open(sys.argv[1]) as fp:
	for line in fp:
		if line.startswith("#CHROM"):
			line = line.rstrip().split()
			for i in range(9, len(line)):
				ind_id[i] = line[i]
				id_N_sites[line[i]] = ""
		elif not line.startswith("#"):
			line = line.rstrip().split()
			variant = line[2]
			for i in range(9, len(line)):
				count = 0
				if line[i] == "0/1":
					count = 1	
				elif line[i] == "1/1":
					count = 2

				if count != 0:
					id_N_sites[ind_id[i]] += variant + "|" + str(count) + ";"

#PMBB_ID	Variants;	Counts;
print("PMBB_ID\tVariant|Count")
for ids in id_N_sites:
	if id_N_sites[ids] != "":
		print(ids + "\t" + id_N_sites[ids].rstrip(";"))

#Next script:
#PMBB_ID	N_variants_carrier	N_genes_carrier	Genes;	Variants_per_gene;	HL	Degree_HL
