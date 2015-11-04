#! /usr/bin/python
'''Converts file in gff to bed format'''
import sys;
import argparse

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts file in gff to bed format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to a file in gff format");
parser.add_argument('-a', '--attributes', nargs = '+', required = False, type = str, help = "Names of gff attributes to be written as additional bed fields in an order they are provided");
args = parser.parse_args();
attrs = args.attributes

for interval in BedTool(args.path):
	if(attrs):
		print "%s\t%d\t%d\t%s\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.stop, interval.name, interval.score, interval.strand, "\t".join([interval.attrs[x] for x in attrs]));
	else:
		print "%s\t%d\t%d\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.stop, interval.name, interval.score, interval.strand);

