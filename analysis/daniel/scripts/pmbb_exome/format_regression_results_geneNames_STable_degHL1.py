import sys

#gene list
gene_list = []
with open(sys.argv[1]) as fp: #include_degHL1/chr22_toModel.txt
	gene_list = fp.readline().split()[26:]

gene_carriers = {}
with open(sys.argv[2]) as fp: #include_degHL1/n_carriers_per_gene_chr22.txt
	fp.readline()
	for line in fp:
		gene, cases, carrier, ctrls = line.rstrip().split()
		gene_carriers[gene] = "\t" + carrier + "\t" + cases

print("Gene\tbeta\tSE\tp\tCarriers\tCase_carriers") #include_degHL1/chr22_linear_notFormatted.txt
count = 0
with open(sys.argv[3]) as fp: 
	for line in fp:
		old_l = line
		line = line.rstrip().split()
		if len(line) > 0:
			if old_l.startswith("d[, i]"):
				p = line[5].strip("<")
				if p == "":
					p = line[6]
				if p != "NA":
					print(gene_list[count] + "\t" + line[2] + "\t" + line[3] +  "\t" + p + "\t" + gene_carriers[gene_list[count]])
				count += 1
