#! /usr/lib/python
'''Changes name for each bed file interval to a new format (example: chr1|+|11873|12227)'''
import argparse
import sys;

from pybedtools import BedTool, Interval



parser = argparse.ArgumentParser(description='Changes name for each bed file interval to a new format (example: chr1:11873-12227(+))');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file");
args = parser.parse_args();

for i in BedTool(args.path):
	if(not i.strand):
		strand = i[3];
	else:
		strand = i.strand
	if(not i.score):
		score = '0'
	else:
		score = i.score
	name = "%s|%s|%d|%d" % (i.chrom, strand, i.start, i.stop)
	j = Interval(i.chrom, i.start, i.stop, name, score, strand)
	sys.stdout.write(str(j))