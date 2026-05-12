import sys

gene = ""
print("Gene\tp")
with open(sys.argv[1]) as fp: 
	for line in fp:
		old_l = line
		line = line.rstrip().split()
		if len(line) > 0:
			if old_l.startswith('[1] "'):
				gene = line[1].strip('"')
			if old_l.startswith("d[, i]"):
				p = line[5].strip("<")
				if p == "":
					p = line[6]
				print(gene + "\t" + p)
