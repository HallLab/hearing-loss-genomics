import sys

hup_pmbb = {}
hup_dob = {}
with open(sys.argv[1]) as fp: #audbase_1.7.21_rj_02192021.csv
	fp.readline()
	for line in fp:
		line = line.rstrip().split(",")

		hup =  line[1]
		pmbb = line[-1]
		hup_pmbb[hup] =  pmbb
		hup_dob[hup] = line[2].strip('"')

with open(sys.argv[2]) as fp: #20210222_AudBase.TXT
	print(fp.readline().rstrip() + ",PMBB_ID")
	for line in fp:
		old_l = line.rstrip()
		line = old_l.split(",")

		hup = line[4].strip('"').lstrip("0")

		day, mon, year = line[5].split("/")
		day = day.lstrip("0")
		mon = mon.lstrip("0")

		date = day + "/" + mon + "/" + year

		if hup in hup_pmbb and hup_dob[hup] == date:
			print(old_l + "," + hup_pmbb[hup])
		else:
			print(old_l + ",NA")
