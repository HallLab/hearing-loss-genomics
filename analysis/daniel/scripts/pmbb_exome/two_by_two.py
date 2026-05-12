import sys

hl = set()
with open(sys.argv[1]) as fp: #gene_list_regions.txt
	for line in fp:
		hl.add(line.rstrip().split()[1])

in_hl = 0
hl_p05 = 0
any_p05 = 0
tot_genes = 0
with open(sys.argv[2]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		p = line[3]
		if p != "NA":
			p = float(line[3])

			if line[0] in hl:
				in_hl += 1
	
				if p < .05:
					hl_p05 += 1

				hl.remove(line[0])

			if p < .05:
				any_p05 += 1

			tot_genes += 1
print(hl)

print(str(hl_p05) + "\t" + str(in_hl - hl_p05))
print(str(any_p05 - hl_p05) + "\t" + str( (tot_genes - in_hl) - (any_p05 - hl_p05)))
