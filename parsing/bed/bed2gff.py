#! /usr/bin/python
'''Converts file in bed to gff format'''
import sys;
import argparse

from pybedtools import BedTool

from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Converts file in bed to gff format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to a file in bed format");
parser.add_argument('--feature', nargs = '?', default='.', type = str, help = "Type of the features in bed file. It corresponds to \'feature\' field in gff format");
args = parser.parse_args();

for interval in BedTool(args.path):
	sys.stdout.write(str(construct_gff_interval(interval.chrom, interval.start, interval.stop, args.feature, score=interval.score, strand=interval.strand, attrs=[('ID', interval.name)])))
	#print "%s\t%d\t%d\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.stop, interval.name, interval.score, interval.strand)