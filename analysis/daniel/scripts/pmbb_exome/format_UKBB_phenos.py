import sys

print("ID\tHL_case\tGender\tAge\tAge_sq\tPC1\tPC2\tPC3\tPC4\tPC5")
with open(sys.argv[1]) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		hl = 0
		if line[2] == "99" or (line[2] == "1" and line[3] == "1"):
			hl = 1
		elif line[2] != line[3]:
			hl = "NA"
		age_sq = "NA"
		try:
			age_sq = str( float(line[4]) * float(line[4]))
		except:
			pass
		print(line[0] + "\t" + str(hl) + "\t" + line[1] + "\t" + line[4] + "\t"  + age_sq + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\t" + line[9])
