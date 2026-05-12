import sys

headers = []
gene_count = {}
with open(sys.argv[1]) as fp:
	headers = fp.readline().rstrip().split()
	for line in fp:
		line = line.rstrip().split()
	
		for i in range(26, len(headers)):
			if headers[i] not in gene_count and line[i] != "NA":
				gene_count[headers[i]] = 0
			if line[i] != "0" and line[i] != "NA" and line[24] != "0":
				gene_count[headers[i]] += 1

print("Gene\tN_carriers_cases")
for g in sorted(gene_count):
	print(g  +"\t" + str(gene_count[g]))
