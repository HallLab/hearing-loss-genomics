import sys 

#parse list of removed cases
#look them up in relatedness file
	#add all their relatives to list to remove
#get original list to remove
#add everything except hl cases
#print the above + hl relatives

#parse list of removed cases
removed = set()
with open(sys.argv[1]) as fp: #removed_cases.txt
	for line in fp:
		removed.add(line.rstrip().split()[0])

#look them up in relatedness file
	#add all their relatives to list to remove
toremove = set()
with open(sys.argv[2]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Exome/IBD/PMBB-Release-2020-2.0_genetic_exome.genome
	fp.readline()
	for line in fp:
		line = line.rstrip().split()

		if line[0] in removed:
			toremove.add(line[1])
		elif line[1] in removed:
			toremove.add(line[0])

#get original list of all unrelated
	#add back removed hl cases
	#remove toremove()
#print the above + hl relatives
tokeep = set()
with open(sys.argv[3]) as fp: #/project/PMBB/PMBB-Release-2020-2.0/Exome/IBD/PMBB-Release-2020-2.0_genetic_exome_3rd_degree_unrelated
	for line in fp:
		tokeep.add(line.rstrip())

for i in removed:
	tokeep.add(i)
for i in toremove:
	if i in tokeep:
		tokeep.remove(i)

for i in sorted(tokeep):
	print(i)

#check if has all hl cases
