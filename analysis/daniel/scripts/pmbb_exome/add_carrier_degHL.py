import sys

pmbb_degHL = {}
with open(sys.argv[1]) as fp: #RGC21_45k_aud_1.csv
	fp.readline()
	for line in fp:
		line = line.rstrip().split(",")
		if line[4].replace('"', '').startswith("PMBB"):
			pmbb_degHL[line[4].replace('"', '')] = line[-2]

cases = set()
with open(sys.argv[2]) as fp: #cases_control.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if str(line[1]) == "1":
			cases.add(line[0])

with open(sys.argv[3]) as fp: #allGenes/HL_rmAudNA/n_carriers_nVariants_rmAudNA.txt
	print("PMBB_ID\tN_HL_genes_carrier\tN_variants_pred_deleterious\tCarrier\tDeg_HL")
	for line in fp:
		if not line.startswith("#") and not line.startswith("PMBB_ID"):
			line = line.rstrip().split()
	
			#carrier or not	
			carrier = 0
			if int(line[1]) > 0:
				carrier = 1

			#degree HL and print
			if line[0] in pmbb_degHL:
				print(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + str(carrier) + "\t" + pmbb_degHL[line[0]])
			elif line[0] not in pmbb_degHL and line[0] not in cases:
				print(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + str(carrier) + "\t0")
			else:
				print(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + str(carrier) + "\tNA")
