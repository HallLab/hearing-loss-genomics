import sys

pmbb_covs = {}
with open(sys.argv[1]) as fp: #rgc21_45k_aud_1.csv
	fp.readline()
	for line in fp:
		line = line.rstrip().split(",")
		birthyear = line[2].strip('"').split("/")[2]
		audyear = line[3].split("-")[0]

		if audyear != "NA" and birthyear != "NA":
			age = float(audyear) - float(birthyear)
			pmbb = line[4].strip('"')
			pmbb_covs[pmbb] = str(age) + "\t" + str(age*age)

print("PMBB_ID\tSex\tAge\tAgeSq\tPC1\tPC2\tPC3\tPC4")
with open(sys.argv[2]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Phenotype/PMBB-Release-2020-2.0_phenotype_covariates.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()

		sex = "NA"
		if line[1] == "Male":
			sex = "1"
		elif line[1] == "Female":
			sex = "0"	

		if line[0] in pmbb_covs:
			print(line[0] + "\t" + sex + "\t" + pmbb_covs[line[0]] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7])
		elif line[3] != "NA":
			print(line[0] + "\t" + sex + "\t" + line[3] + "\t" + str(float(line[3])*float(line[3])) + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7])
