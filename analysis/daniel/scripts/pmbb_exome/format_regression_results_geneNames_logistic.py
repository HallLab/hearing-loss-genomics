import sys

#gene list
gene_list = []
with open(sys.argv[1]) as fp: #allGenes/20PCs/degreeHL/chr22_toModel.txt
	gene_list = fp.readline().split()[27:]

print("Gene\tp")
count = 0
with open(sys.argv[2]) as fp: 
	for line in fp:
		old_l = line
		line = line.rstrip().split()
		if len(line) > 0:
			if old_l.startswith("d[, i]"):
				p = line[5].strip("<")
				if p == "":
					p = line[6]
				print(gene_list[count] + "\t" + p)
				count += 1
