import sys

ind_id = {}
id_N_sites = {}
with open(sys.argv[1]) as fp:
	for line in fp:
		if line.startswith("#CHROM"):
			line = line.rstrip().split()
			for i in range(9, len(line)):
				ind_id[i] = line[i]
				id_N_sites[line[i]] = ""
		elif not line.startswith("#"):
			line = line.rstrip().split()
			for i in range(9, len(line)):
				if line[i] != "0/0" and line[i] != "./.":
					id_N_sites[ind_id[i]] += line[2] + "\t" + line[i] + "\t" 

print("PMBB_ID\tN_sites")
for ids in id_N_sites:
	if id_N_sites[ids] != "":
		print(ids + "\t" + id_N_sites[ids])
