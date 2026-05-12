import sys
import gzip

func_refGene = {"exonic", "exonic;splicing"}
exonicFunc_refGene = {"frameshift substitution", "stopgain"}
with gzip.open(sys.argv[1]) as fp: #annot_genes_full.txt
	print(fp.readline().rstrip())
	for line in fp:
		old_l =  line.rstrip()
		line = old_l.split("\t")

		if (line[87] == "." or line[82] == "." or (float(line[87]) < .001 and float(line[82]) < .001) ) and (line[6] == "splicing" or (line[6] in func_refGene and line[9] in exonicFunc_refGene)) or (line[9] == "nonsynonymous SNV" and line[46] != "."  and float(line[46]) > .60):
			print(old_l)
