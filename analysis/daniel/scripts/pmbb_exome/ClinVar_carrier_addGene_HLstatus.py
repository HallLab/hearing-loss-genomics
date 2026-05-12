import sys


#Next script:
#PMBB_ID        N_variants_carrier      N_genes_carrier Genes;  Variants_per_gene;      HL      Degree_HL

var_gene = {}
with open(sys.argv[1]) as fp: #allGenes/annot_genes_full_funcToInclude_allGenes_ClinVar_onlyHL.txt
	for line in fp:
		line = line.rstrip().split()
		var_gene[line[0]] = line[7]

gene_pattern = {}
with open(sys.argv[2]) as fp: #Hearing\ loss\ genes.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		gene_pattern[line[0]] = line[1]

#HL status
pmbb_cc = {}
with open(sys.argv[3]) as fp: #cases_control.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		pmbb_cc[line[0]] = line[1]

#degree HL, and aud or not
pmbb_deg = {}
with open(sys.argv[4]) as fp: #audbase_feb252021/RGC21_45k_aud_1.csv
	fp.readline()
	for line in fp:
		line = line.rstrip().split(",")
		pmbb = line[4].strip('"')
		if pmbb.startswith("PMBB"):
			pmbb_deg[pmbb] = line[-2] + "\t1"

print("PMBB\tVariant|Count\tGene|Pattern\tHL_case\tDegree_HL\tHasAudiogram")
with open(sys.argv[5]) as fp: #allGenes/ClinVar_carriers/carriers_not_formatted.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		
		variants = line[1].split(";")
		gene = ""
		for v in variants:
			vid = v.split("|")[0]
			if var_gene[vid] in gene_pattern:
				gene += var_gene[vid] + "|" + gene_pattern[var_gene[vid]] + ";"
			else:
				gene += var_gene[vid] + "|NA;"

		if line[0] in pmbb_cc:
			if line[0] in pmbb_deg:
				print(line[0] + "\t" + line[1] + "\t" + gene.rstrip(";") + "\t" + pmbb_cc[line[0]] + "\t" + pmbb_deg[line[0]])
			else:
				print(line[0] + "\t" + line[1] + "\t" + gene.rstrip(";") + "\t" + pmbb_cc[line[0]] + "\tNA\t0")
