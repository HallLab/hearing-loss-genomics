import sys


counts = []
with open(sys.argv[1]) as fp: #e.g. controls_category1.vcf
	for line in fp:
		if not line.startswith("##"):
			if line.startswith("#CHROM"):
				for i in range(0, len(line.rstrip().split()[9:])):
					counts.append(0)
			else:
				line = line.rstrip().split()
				
				for i in range(9, len(line)):
					if "1" in line[i]:
						counts[i-9] += 1

n_var = 0.0
for c in counts:
	if c != 0:
		n_var += 1.0

tot_n = float(len(counts))
perc_with = str(100.0*(n_var / tot_n))
perc_without = str(100.0*((tot_n - n_var) / tot_n))

print("Category\tN_with\tN_without\tN_total\tPerc_with\tPerc_without")
print(sys.argv[1].replace(".vcf", "") + "\t" + str(n_var) + "\t" + str(tot_n - n_var) + "\t" + str(tot_n)  + "\t" + perc_with + "\t" + perc_without)
