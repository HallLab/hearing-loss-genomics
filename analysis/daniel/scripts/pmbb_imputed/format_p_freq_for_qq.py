import sys
import gzip

model = sys.argv[3]

snp = {}
print("SNP\tMAF\tP")
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		snp[line[1]] = line[4]

with gzip.open(sys.argv[2]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()

		if line[7] == "ADD" or line[6] == "ADD":
			p = "NA"
			if model == "linear":
				p = line[11]
			elif model == "binary":
				p = line[12]

			print(line[2] + "\t" + snp[line[2]] + "\t" + p)
		
