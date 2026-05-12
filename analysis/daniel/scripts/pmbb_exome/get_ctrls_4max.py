import sys

non_cases = {}
cases = {}

with open(sys.argv[1]) as fp: #ZNF175/ZNF175_pLOF_Joes.txt
	fp.readline()
	for line in fp:
		cases[line.rstrip()] = {}
		cases[line.rstrip()]["n_ctrls"] = 0

#/project/PMBB/PMBB-Release-2020-2.0/Exome/PCA/PMBB-Release-2020-2.0_genetic_exome_ancestries.txt
with open(sys.argv[2]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] not in cases:
			non_cases[line[0]] =  {}
			non_cases[line[0]]["subpop"] = line[1]
			non_cases[line[0]]["superpop"] = line[2]
		else:
			cases[line[0]]["subpop"] = line[1]
			cases[line[0]]["superpop"] = line[2]

with open(sys.argv[3]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Phenotype/PMBB-Release-2020-2.0_phenotype_covariates.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		
		if line[0] not in cases:
			non_cases[line[0]]["sex"] = line[1]
			non_cases[line[0]]["yob"] = line[2]
		else:
			cases[line[0]]["sex"] = line[1]
			cases[line[0]]["yob"] = line[2]
			
print("Control_ID\tCase_ID\tSubpop_control\tSubpop_case\tSuperpop_control\tSuperpop_case\tBirth_year_control\tBirth_year_case\tSex\tBirth_year_diff")
used = set()
for non_case in non_cases:
	for case in cases:
		s1 = non_cases[non_case]["sex"]
		s2 = cases[case]["sex"]

		y1 = non_cases[non_case]["yob"]
		y2 = cases[case]["yob"]

		sp1 = non_cases[non_case]["subpop"]
		sp2 = cases[case]["subpop"]

		if s1 == s2 and y1 == y2 and sp1 == sp2 and cases[case]["n_ctrls"] < 4 and non_case not in used:
			print(non_case + "\t" + case + "\t" + sp1 + "\t" + sp2 + "\t" + non_cases[non_case]["superpop"] + "\t" + cases[case]["superpop"] + "\t" + y1 + "\t" + y2 + "\t" + s1 + "\t0")
			cases[case]["n_ctrls"] += 1
			used.add(non_case)


for case in cases:
	if cases[case]["n_ctrls"] < 4:

		matches = {}
		for non_case in non_cases:
			s1 = non_cases[non_case]["sex"]
			s2 = cases[case]["sex"]

			y1 = non_cases[non_case]["yob"]
			y2 = cases[case]["yob"]

			sp1 = non_cases[non_case]["superpop"]
			sp2 = cases[case]["superpop"]

			if s1 == s2 and sp1 == sp2 and cases[case]["n_ctrls"] < 4:
				matches[non_case] = abs(float(y1) - float(y2))

		for match in sorted(matches, key=matches.get):
			if cases[case]["n_ctrls"] < 4 and match not in used:
				print(match + "\t" + case + "\t" + non_cases[match]["subpop"] + "\t" + cases[case]["subpop"] + "\t" + non_cases[match]["superpop"]  + "\t" + cases[case]["superpop"] + "\t" + non_cases[match]["yob"] + "\t" + cases[case]["yob"] + "\t" + cases[case]["sex"] + "\t" + str(matches[match]))
				cases[case]["n_ctrls"] += 1
				used.add(match)
