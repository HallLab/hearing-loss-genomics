import sys

id_status = {}
with open(sys.argv[1]) as fp: #cases_control.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		id_status[line[0]] = line[1]

print("PMBB_ID\tSNHL")
with open(sys.argv[2]) as fp: #cases_control_allowOnlyPhecode_formatted.txt
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		if line[0] in id_status and id_status[line[0]] == line[1]:
			print(line[0] + "\t" + line[1])
