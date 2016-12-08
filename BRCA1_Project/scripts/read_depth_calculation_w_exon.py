#!/usr/bin/env python
import sys
#argv[1] is the file with the following format, chr17	41196291	41196295	30	1. The file doesn't have header.
total = 0

with open(sys.argv[1]) as f:
	for line in f:
		array = line.split()
		chrom = array[0]
		start = int(array[1])
		end = int(array[2])
		depth = array[3]
		exon = array[4]
		total = total + (end - start)

		while (start < end):
			print(chrom + "\t" + str(start) + "\t" + depth + "\t" + exon)
			start = start + 1
#print(total)