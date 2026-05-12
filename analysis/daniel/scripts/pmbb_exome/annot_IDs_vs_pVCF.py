import sys

chrPos_ID = {}
with open(sys.argv[1]) as fp: #annot_genes_full_funcToInclude.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		chrPos_ID[line[1] + " " + line[2]] = "NA\tNA\tNA\tNA\tNA\t" +  line[0] + "\t" + line[4] + "\t" + line[5]

with open(sys.argv[2]) as fp: #vcf_SNP_IDs_chr22.txt
	for line in fp:
		old_l = line.rstrip()
		line = old_l.split()
		if (line[0] + " " + line[1]) in chrPos_ID:
			annot = chrPos_ID[line[0] + " " + line[1]].split()
			chrPos_ID[line[0] + " " + line[1]] = old_l + "\t" + annot[5] + "\t" + annot[6] + "\t" + annot[7]


print("chrom\tpos\tVCF_ID\tVCF_ref\tVCF_alt\tannot_ID\tannot_ref\tannot_alt")
for c in chrPos_ID:
	print(chrPos_ID[c])
