import sys

pmbbs = set()
with open(sys.argv[1]) as fp:
	for line in fp:
		if line.startswith("PMBB"):
			pmbbs.add(line.rstrip().split(",")[0])

ncase = 0
nctrl = 0
with open(sys.argv[2]) as fp:
	for line in fp:
		pmbb, pheno = line.rstrip().split()
		if pmbb in pmbbs:
			if pheno == "0":
				nctrl += 1
			elif pheno == "1":
				ncase += 1

print("Cases: " + str(ncase))
print("Controls: " +  str(nctrl))
