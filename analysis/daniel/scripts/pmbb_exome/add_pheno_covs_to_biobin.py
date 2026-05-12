import sys

pmbb_covs = {}
with open(sys.argv[1]) as fp: #covs_withAnc_onlyEUR-AFR_rmAncColumn_20PCs.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
			
		toprint = line[0]
		for i in range(1, len(line)):
			toprint += "\t" + line[i]
		
		pmbb_covs[line[0]] = toprint

with open(sys.argv[2]) as fp: #cases_control_degHL.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] in pmbb_covs:
			pmbb_covs[line[0]] += "\t" + line[1] + "\t" + line[2]

header = "PMBB_ID\tPC1\tPC2\tPC3\tPC4\tSex\tAge\tAgeSq\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\tPC11\tPC12\tPC13\tPC14\tPC15\tPC16\tPC17\tPC18\tPC19\tPC20\tDegHL\tDegHL_rand"

with open(sys.argv[3]) as fp: #allGenes/20PCs/degreeHL/chr22-bins.csv
	for line in fp:
		if line.startswith("ID") or line.startswith("PMBB"):
			line = line.rstrip().split(",")
			toprint = ""
			for i in range(2, len(line)):
				toprint += "\t" + line[i]

			#first line
			if line[0].startswith("ID"):
				print(header + toprint)
			#others
			elif line[0] in pmbb_covs:
				print(pmbb_covs[line[0]] + toprint)
	
			

		
