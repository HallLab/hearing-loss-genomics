import sys

pmbb_pcs = {}
with open(sys.argv[1]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Genotype/PCA/PMBB-Release-2020-2.0_genetic_genotype.eigenvec
	fp.readline()
	for line in fp:
		old_l = line.rstrip()
		pmbb = old_l.split()[0]

		pmbb_pcs[pmbb] = old_l

print("FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\tPC11\tPC12\tPC13\tPC14\tPC15\tPC16\tPC17\tPC18\tPC19\tPC20\tGender\tAge\tAgeSq")
with open(sys.argv[2]) as fp: #../PMBB_Exome/include_degHL1/chr22_toModel.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] in pmbb_pcs:
			print(line[0] + "\t" + pmbb_pcs[line[0]] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7])
