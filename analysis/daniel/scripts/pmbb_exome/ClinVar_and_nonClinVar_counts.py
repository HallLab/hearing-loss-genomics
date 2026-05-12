import sys

cv = set()
with open(sys.argv[1]) as fp: #allGenes/annot_genes_full_funcToInclude_allGenes_ClinVar_onlyHL.extract
	for line in fp:
		cv.add(line.rstrip())
		

ind_id = {} #index mapping to ID
id_N_sites = {} #ID mapping to N sites
with open(sys.argv[2]) as fp:
	for line in fp:
		if line.startswith("#CHROM"):
			line = line.rstrip().split()
			for i in range(9, len(line)):
				ind_id[i] = line[i]
				id_N_sites[line[i]] = {}
				id_N_sites[line[i]]["clinvar"] = 0
				id_N_sites[line[i]]["non_clinvar"] = 0
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
					if variant in cv:
						id_N_sites[ind_id[i]]["clinvar"] += count
					else:
						id_N_sites[ind_id[i]]["non_clinvar"] += count

print("PMBB_ID\tClinVar_count\tNon_ClinVar_count")
for ids in id_N_sites:
	print(ids + "\t" + str(id_N_sites[ids]["clinvar"]) + "\t" + str(id_N_sites[ids]["non_clinvar"]))
