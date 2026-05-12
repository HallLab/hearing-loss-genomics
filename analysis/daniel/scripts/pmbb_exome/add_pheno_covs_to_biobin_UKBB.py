import sys

ukbb_covs = {}
with open(sys.argv[1]) as fp: #covs.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
			
		toprint = line[0]
		for i in range(1, len(line)):
			toprint += "\t" + line[i]
		
		ukbb_covs[line[0]] = toprint


header = "ID_ID\tGender\tAge\tAge_sq\tPC1\tPC2\tPC3\tPC4\tPC5\tHL"


with open(sys.argv[2]) as fp: #UKBB_analyses/pheno_hearing_aid.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		
		if line[0] in ukbb_covs:
			ukbb_covs[line[0]] += "\t" + line[1]


with open(sys.argv[3]) as fp: #results_hearing_aid/chr22-bins.csv
	for line in fp:
		if line.startswith("ID") or line[0].isdigit():
			line = line.rstrip().split(",")
			toprint = ""
			for i in range(2, len(line)):
				toprint += "\t" + line[i]

			#first line
			if line[0].startswith("ID"):
				print(header + toprint)
			#others
			elif line[0] in ukbb_covs:
				print(ukbb_covs[line[0]] + toprint)
	
			

		
