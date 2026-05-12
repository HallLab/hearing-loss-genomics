import sys

snp_annots = {}
with open(sys.argv[1]) as fp: #annot_toInclude_TCOF1.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split("\t")

		snp_annots[line[0]] = line[6] + "\t" + line[9] + "\t" + line[51] + "\t" + line[99]

pmbb_hl = {}
with open(sys.argv[2]) as fp: #include_degHL1/cases_control_degHL1.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		pmbb_hl[line[0]] = line[2]

print("PMBB\tDegHL\tVariant\tVariant_count\tFunc.refGene\tExonicFunc.refGene\tREVEL\tClinVar")
with open(sys.argv[3]) as fp: #TCOF1_carriers_studyIndvs.txt
	for line in fp:
		line = line.rstrip().split()

		print(line[0].split("_")[0] + "\t" + pmbb_hl[line[0].split("_")[0]] + "\t" + line[1] + "\t" + line[2] + "\t" + snp_annots[line[1]])
