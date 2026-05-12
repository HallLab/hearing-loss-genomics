import sys
import gzip

model = sys.argv[3]

snp = {}
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		snp[line[1]] = line[4]

if model == "linear":
	print("CHROM\tPOS\tSNP\tMAF\tBETA\tSE\tP")
elif model == "binary":
	
with gzip.open(sys.argv[2]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()

		if line[7] == "ADD" or line[6] == "ADD":
			p = "NA"
			beta = "NA"
			se = "NA"
			chrom = line[0]
			pos = line[1]

			if model == "linear":
				p = line[11]
				beta = line[8]
				se = line[9]
			elif model == "binary":
				p = line[12]

			print(chrom +"\t" + pos + "\t" +line[2] + "\t" + snp[line[2]] + "\t" +beta + "\t" + se +"\t" + p)
