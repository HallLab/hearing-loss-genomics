import sys


iid_SNPs = {}
header = []
with open(sys.argv[1]) as fp:
	for line in fp:
		if not line.startswith("##"):
			if line.startswith("#"):
				header = line.rstrip().split()
			else:
				line = line.rstrip().split()
				for i in range(0, len(line)):
					if line[i] == "1/1":
						if header[i] not in iid_SNPs:
							iid_SNPs[header[i]] = ""
						iid_SNPs[header[i]] += "\t" + line[2]

cv = set()
with open(sys.argv[2]) as fp: #clinvar
	fp.readline()
	for line in fp:
		cv.add(line.rstrip().split()[0])

for iid in sorted(iid_SNPs):
	if iid_SNPs[iid] in cv:
		print(iid + "\t" + iid_SNPs[iid] + "\tYes")
	else:
		print(iid + "\t" + iid_SNPs[iid] + "\tNo")
