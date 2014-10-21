#! /usr/lib/python
'''assigns genomic to the chimeric reads which were mapped to non-genomic reference(transcriptome, 3'utr, etc.)'''
import argparse;
import sys;

def reassign(a):
	chrom, strand, start, stop = a[0].split("|")
	start = int(start)
	stop = int(stop)
	
	if(strand == '+'):
		stop = start + int(a[2])
		start = start + int(a[1])		
	else:
		start = stop - int(a[2])
		stop = stop - int(a[1])
	
	return chrom, str(start), str(stop), a[3], a[4], strand

for l in sys.stdin:
	a = l.strip().split("\t");
	a[:6] = reassign(a[:6])
	a[6:12] = reassign(a[6:12])
	print "\t".join(a);

#parser = argparse.ArgumentParser(description='produce chimeras from merged and already filtered sam file');
#parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to file with chimeras");
#args = parser.parse_args();