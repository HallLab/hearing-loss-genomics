import sys

case_control = {}
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		case_control[line[0]] = line[1]


headers = []
gene_count = {}
with open(sys.argv[2]) as fp:
	headers = fp.readline().rstrip().split(",")
	for line in fp:
		line = line.rstrip().split(",")

		if "PMBB" in line[0]:	
			for i in range(2, len(headers)):
				#init
				if headers[i] not in gene_count:
					gene_count[headers[i]] = {}
					gene_count[headers[i]]["cases"] = 0
					gene_count[headers[i]]["controls"] = 0
					gene_count[headers[i]]["total"] = 0

				#add if carrier
				line[i] = int(line[i])
				if line[i] > 0:
					if case_control[line[0]] == "0":
						gene_count[headers[i]]["controls"] += 1
					elif case_control[line[0]] == "1":
						gene_count[headers[i]]["cases"] += 1
					gene_count[headers[i]]["total"] += 1
					

print("Gene\tN_carriers_cases\tN_carriers_total\tN_carriers_controls")
for g in sorted(gene_count):
	print(g  +"\t" + str(gene_count[g]["cases"]) + "\t" + str(gene_count[g]["total"]) + "\t" + str(gene_count[g]["controls"]) )
