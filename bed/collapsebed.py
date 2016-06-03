#! /usr/bin/python
'''Collapses intervals with identical coordinates. Provided gff/bed file has to be sorted'''


import sys;
import argparse

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Collapses intervals with identical coordinates. Provided gff/bed file has to be sorted');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic intervals. Provided gff/bed file has to be sorted");
args = parser.parse_args();


curcoords = ();
count = 0;
for interval in BedTool(args.path):
	coords = (interval.chrom, interval.strand, interval.start, interval.end)
	if(coords != curcoords):
		count += 1
		interval.name = "collapsed_%d" % count
		sys.stdout.write(str(interval))
		curcoords = coords