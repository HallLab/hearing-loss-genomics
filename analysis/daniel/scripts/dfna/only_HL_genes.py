import sys

genes = set()
with open(sys.argv[1]) as fp: #postlingual_HL_genes.txt
	fp.readline()
	for line in fp:
		genes.add(line.rstrip().split()[0])

with open(sys.argv[2]) as fp: #annot_only_exonic_splicing.txt
	print(fp.readline().rstrip())
	for line in fp:
		old_l = line.rstrip()
		line = old_l.split()

		if line[6] in genes:
			print(old_l)
