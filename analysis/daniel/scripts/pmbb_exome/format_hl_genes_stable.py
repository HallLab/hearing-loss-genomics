import sys

genes = {}
with open(sys.argv[1]) as fp: #Hearing_loss_genes.txt
	for line in fp:
		line = line.rstrip().split()
		genes[line[0]] = line[1]

print("Gene\tInheritance_pattern\tChromosome_hg38")
with open(sys.argv[2]) as fp: #gene_list_regions.txt
	for line in fp:
		line = line.rstrip().split()

		if line[0] != "X":
			if line[1] in genes:
				print(line[1] + "\t" + genes[line[1]] + "\t" + line[0])
			else:
				print(line[1] + "\tNA\t" + line[0])
