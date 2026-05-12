import sys

hup_pmbb = {}
with open(sys.argv[1]) as fp: #audbase_1.7.21_rj_02192021.csv
	fp.readline()
	for line in fp:
		line = line.rstrip().split(",")

		hup =  line[1]
		pmbb = line[-1]
		hup_pmbb[hup] =  pmbb

with open(sys.argv[2]) as fp: #20210222_AudBase.TXT
	print(fp.readline().rstrip() + ",PMBB_ID")
	for line in fp:
		old_l = line.rstrip()
		line = old_l.split(",")

		hup = line[4].strip('"').lstrip("0")

		if hup in hup_pmbb:
			print(old_l + "," + hup_pmbb[hup])
		else:
			print(old_l + ",NA")
