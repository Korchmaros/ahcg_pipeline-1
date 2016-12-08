import sys

with open (sys.argv[1]) as f_nm:
	nm = f_nm.readlines()
	with open (sys.argv[2]) as f_final_bed:
		final_bed = f_final_bed.readlines()

for bed_range in final_bed:
	for nm_range in nm:
		nm_range_left = int(nm_range.split("\t")[1])
		nm_range_right = int(nm_range.split("\t")[2])
		bed_range_left = int(bed_range.split("\t")[1])
		bed_range_right = int(bed_range.split("\t")[2])
		exon_number = nm_range.split("\t")[-1]
		if nm_range_left <= bed_range_left and bed_range_right <= nm_range_right:
			print(str(bed_range).rstrip() + "\t" + str(exon_number).rstrip())