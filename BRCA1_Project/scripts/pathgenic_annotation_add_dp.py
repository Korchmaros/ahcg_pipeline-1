#! /usr/local/env python3.4
#argv[1] is the vcf file 
#argv[2] is the csv annotation file download from brca database
import csv
import vcf
import sys

vcff  =  open(sys.argv[1])
csvf = sys.argv[2]

csv_dict = {lines[-7]  : lines[14] for count, lines in enumerate(csv.reader(open(csvf), delimiter=',')) if count != 0}

for lines in vcf.Reader(vcff):
	if '{0}:{1}:{2}>{3}'.format(lines.CHROM, lines.POS, lines.REF, str(lines.ALT[0])) in csv_dict:
		print( lines,csv_dict[ '{0}:{1}:{2}>{3}'.format(lines.CHROM, lines.POS, lines.REF, str(lines.ALT[0]))], lines.INFO["DP"])