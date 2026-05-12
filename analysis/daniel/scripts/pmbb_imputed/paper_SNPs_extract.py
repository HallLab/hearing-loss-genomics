import sys

all_cp = set()
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split("\t")
		chrom_pos = line[4] + " " + line[5]
		all_cp.add(chrom_pos)

with open(sys.argv[2]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		chrom_pos = line[0] + " " + line[1]
		if chrom_pos in all_cp:
			print(line[2])
