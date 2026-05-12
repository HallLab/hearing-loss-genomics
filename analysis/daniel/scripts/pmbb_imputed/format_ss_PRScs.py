import sys

chrpos_id = {}
with open(sys.argv[1]) as fp: #hm3_allchr.bim
	for line in fp:
		line = line.rstrip().split()
		chrpos_id[line[0] + " " + line[3]] = line[1]

print("SNP\tA1\tA2\tBETA\tP")
with open(sys.argv[2]) as fp: #GWAS_hm3_hg19_formatted.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		
		chrpos = line[0] + " " + line[1]
		
		if chrpos in chrpos_id:
			print(chrpos_id[chrpos] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[-1])
