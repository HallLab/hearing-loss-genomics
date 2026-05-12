import sys

pmbb_mrn = {}
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split(",")
		pmbb_mrn[line[-1].replace('"', '')] = line[4].replace('"', '')

with open(sys.argv[2]) as fp:
	print("MRN\tPMBB")
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		pmbb = line[0].split("_")[0]

		if pmbb in pmbb_mrn:
			print(pmbb_mrn[pmbb] + "\t" + pmbb)
		else:
			print("NA\t" + pmbb)
