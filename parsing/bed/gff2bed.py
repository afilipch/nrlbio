#! /usr/bin/python
'''Converts file in gff to bed format'''
import sys;
import argparse

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts file in gff to bed format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to a file in gff format");
args = parser.parse_args();

for interval in BedTool(args.path):
	print "%s\t%d\t%d\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.stop, interval.name, interval.score, interval.strand)

