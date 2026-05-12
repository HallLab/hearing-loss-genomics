import sys

pmbb_anc = {}
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		pmbb_anc[line[0]] = line[2]

with open(sys.argv[2]) as fp:
	print(fp.readline().rstrip() + "\tAncestry")
	for line in fp:
		old_l = line.rstrip()
		line = old_l.split()
		print(old_l + "\t" + pmbb_anc[line[0]])
