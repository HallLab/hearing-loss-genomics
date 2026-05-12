import sys

indir = sys.argv[1]

print("Chromosome_hg38\tGene\tLogistic_regression_beta_p")
for i in range(1, 23):
	with open(indir + "/chr" + str(i) + "_results_t.txt") as fp:
		fp.readline()
		for line in fp:
			line = line.rstrip().split("\t")
			if line[1] != "NA" and line[8] != "NA" and line[8] != "nan":
				print(str(i) + "\t" + line[0] + "\t" + line[8])
