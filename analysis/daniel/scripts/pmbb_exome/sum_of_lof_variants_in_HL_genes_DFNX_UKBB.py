import sys

prefix = sys.argv[2] #HL_needAud_merged_maf.01_noRels_keepHLcases_allGenes_chr
#postfix = sys.argv[3] #binary or linear

#list of hl genes
hl_genes = set()
with open(sys.argv[1]) as fp:#gene_list_regions.txt
	for line in fp:
		hl_genes.add(line.rstrip())

found = set()
#counts of:
#-hl genes with >=1 lof variant
#-sum of all lof variants in hl genes
genes = ""
pmbb_counts = {}
for i in range(1, 23):
#for i in range(1, 2):
	with open(prefix + str(i) + "-bins.csv") as fp:
		genes = fp.readline().rstrip().split(",")
		for line in fp:
			line = line.rstrip().split(",") 

			pmbb = line[0]

			for i in range(1, len(genes)):
				if genes[i] in hl_genes and pmbb[0].isdigit():
					found.add(genes[i])
					if pmbb not in pmbb_counts:
						pmbb_counts[pmbb] = {}
						pmbb_counts[pmbb]["n_genes"] = 0
						pmbb_counts[pmbb]["n_variants"] = 0
					if int(line[i]) > 0:
						pmbb_counts[pmbb]["n_genes"] += 1
						pmbb_counts[pmbb]["n_variants"] += int(line[i])

sym_diff = hl_genes^found
for i in sorted(sym_diff):
	print("#" + i)

#print out
print("PMBB_ID\tN_HL_genes_carrier\tN_variants_pred_deleterious")
for pmbb in sorted(pmbb_counts):
	print(pmbb + "\t" + str(pmbb_counts[pmbb]["n_genes"]) + "\t" + str(pmbb_counts[pmbb]["n_variants"]))
