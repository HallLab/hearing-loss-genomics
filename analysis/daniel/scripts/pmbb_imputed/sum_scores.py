import sys

scores = sys.argv[1]

id_score = {}
with open(scores) as fp:
	for line in fp:
		line = line.rstrip().split()

		if line[0] != "FID":
	
			iid = line[1]
			score = float(line[-1])

			if iid not in id_score:
				id_score[iid] = score
			else:
				id_score[iid] += score

print("IID\tPRS")
for iid in sorted(id_score):
	print(iid + "\t" + str(id_score[iid]))
