import sys

headers = []
gene_count = {}
with open(sys.argv[1]) as fp:
	headers = fp.readline().rstrip().split()
	for line in fp:
		line = line.rstrip().split()
	
		for i in range(27, len(headers)):
			if line[i] != "0" and line[i] != "NA":
				if headers[i] not in gene_count:
					gene_count[headers[i]] = 0
				gene_count[headers[i]] += 1

print("Gene\tN_carriers")
for g in sorted(gene_count):
	print(g  +"\t" + str(gene_count[g]))
