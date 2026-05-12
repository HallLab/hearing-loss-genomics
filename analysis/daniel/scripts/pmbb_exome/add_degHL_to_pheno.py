import sys

pmbb_degHL = {}
exclude = set()
with open(sys.argv[1]) as fp: #RGC21_45k_aud_1.csv
	fp.readline()
	for line in fp:
		line = line.rstrip().split(",")
		iid = line[4].replace('"', '')
		if iid.startswith("PMBB"):
			if line[-1] == "NA":
				exclude.add(iid)
			else:
				pmbb_degHL[iid] = line[-2]

cases = set()
print("PMBB_ID\tDeg_HL")
with open(sys.argv[2]) as fp: #cases_control.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] in exclude:
			print(line[0] + "\tNA")
		else:
			if line[0] in pmbb_degHL:
				print(line[0] + "\t" + pmbb_degHL[line[0]])
			else:
				print(line[0] + "\t0")
