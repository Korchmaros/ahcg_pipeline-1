import sys
#argv[1] is the file which contains two columns. The first column contains name of gene. The 2nd column contains the NM_. The file has header.
#argv[2] is human reference genome file, which you could extract the genome coordinates by given the NM_ number.

def header():
	print("chr"+"\t"+"chrom_start"+"\t"+"chrom_end"+"\t"+"strand"+"\t"+"name")

with open(sys.argv[1]) as fin:
		next(fin)
		match = []
		for line_fin in fin:
			with open(sys.argv[2]) as fin_ref:
				for line_ref in fin_ref:
					if line_fin.split()[1].rstrip() in line_ref:
						match.append(line_ref)				
header()
for line in match:
	chrom_start = line.split()[9].split(",")
	chrom_end = line.split()[10].split(",")
	for x in range(0,len(chrom_start)-1):
		result = line.split()[2]+"\t"+chrom_start[x]+"\t"+chrom_end[x]+"\t"+line.split()[3]+"\t"+line.split()[1]
		print(result)