#! /usr/lib/python
'''assigns genomic to the chimeric reads which were mapped to non-genomic reference(transcriptome, 3'utr, etc.)'''
import argparse;
import sys;

parser = argparse.ArgumentParser(description='assigns genomic coordinates to the regions which are on non-genomic reference(transcriptome, 3\'utr, etc.)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed/gff file");
args = parser.parse_args();

def reassign(a):
	chrom, strand, start, stop = a[0].split("|")[:4]
	start = int(start)
	stop = int(stop)
	
	if(strand == '+'):
		stop = start + int(a[2])
		start = start + int(a[1])		
	else:
		start = stop - int(a[2])
		stop = stop - int(a[1])
	
	return chrom, str(start), str(stop), a[3], a[4], strand

with open(args.path) as f:
	for l in f:
		a = l.strip().split("\t");
		a[:6] = reassign(a[:6]);
		print "\t".join(a);

