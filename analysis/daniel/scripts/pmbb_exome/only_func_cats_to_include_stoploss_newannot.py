import sys

func_refGene = {"exonic", "exonic;splicing"}
exonicFunc_refGene = {"frameshift substitution", "stopgain", "stoploss"}
with open(sys.argv[1]) as fp: #annot_genes_full.txt
	print(fp.readline().rstrip())
	for line in fp:
		old_l =  line.rstrip()
		line = old_l.split("\t")

		if (line[87] == "." or line[92] == "." or (float(line[87]) <= .001 and float(line[92]) <= .001) ) and (line[6] == "splicing" or (line[6] in func_refGene and line[9] in exonicFunc_refGene)) or (line[9] == "nonsynonymous SNV" and line[51] != "."  and float(line[51]) > .60):
			print(old_l)
