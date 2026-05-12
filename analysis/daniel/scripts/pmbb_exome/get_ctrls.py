import sys

ids = {}
with open(sys.argv[1]) as fp: #ZNF175/ZNF175_pLOF_Joes.txt
	fp.readline()
	for line in fp:
		ids[line.rstrip()] = {}
		ids[line.rstrip()]["n_ctrls"] = 0
		ids[line.rstrip()]["subpop"] = "NA"
		ids[line.rstrip()]["superpop"] = "NA"

with open(sys.argv[2]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Exome/PCA/PMBB-Release-2020-2.0_genetic_exome_ancestries.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] in ids:
			ids[line[0]]["subpop"] = line[1]
			ids[line[0]]["superpop"] = line[2]

with open(sys.argv[3]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Phenotype/PMBB-Release-2020-2.0_phenotype_covariates.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] in ids:
			ids[line[0]]["sex"] = line[1]
			ids[line[0]]["yob"] = int(line[2])
			ids[line[0]]["ctrls"] = set()

with open(sys.argv[3]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Phenotype/PMBB-Release-2020-2.0_phenotype_covariates.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		for case in ids:
			if line[0] != case and "sex" in ids[case] and ids[case]["sex"] == line[1] and line[2] != "NA" and (ids[case]["yob"] == int(line[2]) or ids[case]["yob"] == int(line[2]) + 1 or ids[case]["yob"] == int(line[2]) - 1) :
				ids[case]["ctrls"].add(line[0])

#print(ids)

print("Control\tCase\tAncestry\tSex\tBirth_year")
with open(sys.argv[2]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Exome/PCA/PMBB-Release-2020-2.0_genetic_exome_ancestries.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		for case in ids:
			if line[0] in ids[case]["ctrls"] and ids[case]["subpop"] == line[1]:
				print(line[0] + "\t" + case + "\t" + ids[case]["anc"] + "\t" + ids[case]["sex"] + "\t" + str(ids[case]["yob"]))
				ids[case]["n_ctrls"] += 1
