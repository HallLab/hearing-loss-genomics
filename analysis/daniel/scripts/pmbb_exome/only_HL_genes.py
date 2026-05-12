import sys

genes = set()
with open(sys.argv[1]) as fp: #justGenes_formatted.txt
	for line in fp:
		genes.add(line.rstrip())

with open(sys.argv[2]) as fp: #/project/pmbb_all/PMBB-Release-2020-2.0/Exome/Variant_annotations/PMBB-Release-2020-2.0_genetic_exome_variant-annotations.txt
	print(fp.readline().rstrip())
	for line in fp:
		old_l = line.rstrip()
		line = old_l.split()

		if line[7] in genes:
			print(old_l)
