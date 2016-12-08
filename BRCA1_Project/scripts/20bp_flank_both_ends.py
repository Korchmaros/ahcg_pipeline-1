import sys
#argv[1] is the BED file with the following format "chr	chrom_start	chrom_end	strand	name". The file has header.
with open (sys.argv[1]) as fin:
	next(fin)
	for line in fin:
		newline = str(line.split()[0] + "\t" + str(int(line.split()[1]) - 20) + "\t" + str(int(line.split()[2]) + 20) + "\t" + line.split()[3] + "\t" + line.split()[4])
		print(newline)