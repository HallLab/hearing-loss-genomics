import sys

id_hl = {}
with open(sys.argv[1]) as fp: #UKBB_analyses/pheno.txt
	fp.readline()
	for line in fp:
		ID, hl = line.rstrip().split()
		id_hl[ID] = hl

with open(sys.argv[2]) as fp: #UKBB_analyses/hearing_aid_user.txt
	fp.readline()
	for line in fp:
		ID, hl = line.rstrip().split()
		ID = ID + "_" + ID
		
		if ID in id_hl and hl == id_hl[ID]:
			print(ID + "\t" + hl)
	
