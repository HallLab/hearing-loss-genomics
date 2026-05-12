import sys

#case if
	#BL_SNHL is true
#control if 
	#BL_SNHL is false OR
	#phecode is FALSE
#NA if
	#389 is true and missing audiogram
	#what if BL_SNHL is true and phecode is false
		#I think true still


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

print("PMBB_ID\tSNHL")
for p in sorted(pmbb_snhl):
	print(p + "\t" + pmbb_snhl[p])
