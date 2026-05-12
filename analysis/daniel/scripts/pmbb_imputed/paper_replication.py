import sys

all_cp = set()
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split("\t")
		chrom_pos = line[4] + " " + line[5]
		all_cp.add(chrom_pos)

print("Chrom\tPos\tMAF\tBeta\tSE\tP")
with open(sys.argv[2]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		chrom_pos = line[0] + " " + line[1]
		if chrom_pos in all_cp:
			print(line[0] +"\t" + line[1] + "\t" + line[3] + "\t" + line[4] + "\t" +line[5] +"\t" + line[6])
