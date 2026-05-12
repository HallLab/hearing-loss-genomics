import sys

#case if
	#phecode == true and (BL_SNHL == true or BL_SNHL == missing)
#control if 
	#phecode == false
#NA if
	#audiogram != missing and 389 != audiogram


pmbb_snhl = {}
with open(sys.argv[1]) as fp: #audbase_feb252021/RGC21_45k_aud_1.csv
	fp.readline()
	for line in fp:
		line = line.rstrip().split(",")

		pmbb = line[4].strip('"')

		if pmbb != "NA":		
			snhl = line[-1].strip('"')
			pmbb_snhl[pmbb] = snhl


with open(sys.argv[2]) as fp: #phecode_hl.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] not in pmbb_snhl:
			if line[1] == "FALSE":
				pmbb_snhl[line[0]] = "0"
			elif line[1] == "TRUE":
				pmbb_snhl[line[0]] = "1"
		else:
			if (pmbb_snhl[line[0]] == "1" and line[1] == "FALSE") or (pmbb_snhl[line[0]] == "0" and line[1] == "TRUE"):
				pmbb_snhl[line[0]] = "NA\tPhecode_Audiogram_Different"
			elif pmbb_snhl[line[0]] == "NA":
				if line[1] == "TRUE":
					pmbb_snhl[line[0]] = "1\tNA" 
				elif line[1] == "FALSE":
					pmbb_snhl[line[0]] = "0\tNA" 
				else:
					pmbb_snhl[line[0]] = "NA\tNA" 

print("PMBB_ID\tSNHL\tReasonIfApplicable")
for p in sorted(pmbb_snhl):
	print(p + "\t" + pmbb_snhl[p])
