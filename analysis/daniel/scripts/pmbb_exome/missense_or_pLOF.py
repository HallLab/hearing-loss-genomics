import sys

snp_func = {}
annot_ZNF = {}
with open(sys.argv[1]) as fp: #ZNF175/ZNF175_annot_genes_full_funcToInclude.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split("\t")

		snp = line[0]
		func = line[9]
		snp_func[snp] = func

with open(sys.argv[2]) as fp: #ZNF175/ZNF175_matched_snp_IDs_annot_pVCF.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()

		snpZNF = line[2]
		snpAnnot = line[5]

		annot_ZNF[snpZNF] = snpAnnot

with open(sys.argv[3]) as fp: #ZNF175/ZNF175_carriers_allIndvs.txt
	print(fp.readline().rstrip() + "\tcategory")
	for line in fp:
		line = line.rstrip().split()
		if snp_func[annot_ZNF[line[1]]] == "nonsynonymous SNV":
			print(line[0] + "\t" + line[1] + "\tmissense")
		else:
			print(line[0] + "\t" + line[1] + "\tpLOF")

