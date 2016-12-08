import sys

with open (sys.argv[1]) as fin:
		nm_array = []
		for line in fin:
			nm_name  = line.split("\t")[4].rstrip()
			nm_array.append(nm_name)
			nm_occurence = nm_array.count(nm_name)
			print(str(line).rstrip() + "\t" + str(nm_occurence))