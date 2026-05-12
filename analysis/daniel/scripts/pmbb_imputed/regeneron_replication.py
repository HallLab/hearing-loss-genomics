import sys

snps = {}
with open(sys.argv[1]) as fp: #media-1\ \(2\).txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split("\t")
		
		snp = line[1].rstrip().split(":")
		chrpos = snp[0] + "\t" + snp[1]

		OR = line[8].strip('"').split()[0]

		af = line[7]

		direction = "NA"
		if float(OR) > 1.0:
			direction = "+"
		elif float(OR) <= 1.0:
			direction = "-"
		
		p = line[9]
		alldir = line[10]

		snps[chrpos] = af + "\t" + OR + "\t" + direction + "\t" + p + "\t" + alldir


print("CHROM\tPOS\tID\tPMBB_MAF\tPMBB_BETA\tPMBB_SE\tPMBB_P\tDirection_match\tRegen_AF\tRegen_OR\tRegen_direction\tRegen_P\tRegen_meta_alldir")
with open(sys.argv[2]) as fp:
	for line in fp:
		if not line.startswith("CHROM"):
			old_l = line.rstrip()
			line = old_l.split()

			chrpos = line[0] +"\t" + line[1]		
	
			if chrpos in snps:
	
				direction = snps[chrpos].split()[2]
		
				match = "No"
				if (float(line[4]) > 0 and direction == "+") or (float(line[4]) < 0 and direction == "-"):
					match = "Yes"
			
				print(old_l + "\t" + match + "\t"+ snps[chrpos])
