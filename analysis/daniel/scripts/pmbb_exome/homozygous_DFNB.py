import sys

dfnb_indices = {}
with open(sys.argv[1]) as fp:
	for line in fp:
		dfnb_indices[line.rstrip()] = ""

print("PMBB_ID\tDFNB_gene\tCount\tColumn")
with open(sys.argv[2]) as fp: #e.g., include_degHL1/linear_chr22-bins_tabs.txt
	#get index of DFNB gene in header
	header = fp.readline().rstrip().split()
	for i in range(0, len(header)):
		if header[i] in dfnb_indices:
			dfnb_indices[header[i]] = i
		
	
	for line in fp:
		line = line.rstrip().split()

		for i in range(1, len(line)):
			if float(line[i]) >= 2 and header[i] in dfnb_indices:
				print(line[0] + "\t" + header[i] + "\t" + line[i] + "\t" + str(i))
