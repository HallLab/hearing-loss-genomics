import sys
import gzip

chrpos_37to38 = {}
with open(sys.argv[1]) as fp: #hapmap3.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		chrpos_37 = line[0] + " " + line[1]
		chrpos_38 = line[0] + " " + line[-1]
		chrpos_37to38[chrpos_38] = chrpos_37 + " " + line[4]

chrpos_snpID = {}
with open(sys.argv[2]) as fp: #all.bim
	for line in fp:
		line = line.rstrip().split()
		chrpos = line[0] + " " + line[3]
		chrom, pos, rsid = chrpos_37to38[chrpos].split()
		print(line[0] +"\t" + rsid + "\t0\t" + pos + "\t" + line[4] + "\t" + line[5])

