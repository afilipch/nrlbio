#! /usr/lib/python
'''Changes name for each bed file interval to a new format (example: chr1|+|11873|12227)'''
import argparse
import sys;

from pybedtools import BedTool



parser = argparse.ArgumentParser(description='Changes name for each bed file interval to a new format (example: chr1:11873-12227(+))');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file");
args = parser.parse_args();

for i in BedTool(args.path):
	i.name = "%s|%s|%d|%d" % (i.chrom, i.strand, i.start, i.stop)
	sys.stdout.write(str(i))