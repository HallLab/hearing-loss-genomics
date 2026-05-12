import sys

id_snhl = {}
with open(sys.argv[1]) as fp: #cases_control.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		id_snhl[line[0]] = line[1]

with open(sys.argv[2]) as fp:
	print(fp.readline().rstrip() + "\tSNHL")
	for line in fp:
		old_l = line.rstrip()
		line = old_l.split()
		if line[1] in id_snhl:
			print(old_l + "\t" + id_snhl[line[1]])
		else:
			print(old_l + "\tNA")
