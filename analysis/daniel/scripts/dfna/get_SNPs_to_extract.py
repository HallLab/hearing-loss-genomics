import sys

chrposallele = set()
with open(sys.argv[1]) as fp: #e.g. category1.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()

		pos = line[1] #Start, not End, for indels
		chrom = line[0]
		a1 = line[3]
		a2 = line[4]

		snp = pos + " " + chrom + " " + a1 + " " + a2
		chrposallele.add(snp)


with open(sys.argv[2]) as fp: #ukb_50k_fe_exome.bim
	for line in fp:
		line = line.rstrip().split()

		pos = line[3] #Start, not End, for indels
		chrom = line[0]
		a1 = line[4]
		a2 = line[5]

		snp1 = pos + " " + chrom + " " + a1 + " " + a2
		snp2 = pos + " " + chrom + " " + a2 + " " + a1

		if snp1 in chrposallele or snp2 in chrposallele:
			print(line[1])
