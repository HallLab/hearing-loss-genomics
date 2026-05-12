import sys

with open(sys.argv[1]) as fp:
	for line in fp:
		if not line.startswith("#"):
			if "R2=" in line:
				r2 = line.rstrip().split("R2=")[1].split(";")[0]
	
				if float(r2) < .30:
					print(line.rstrip().split()[2])
