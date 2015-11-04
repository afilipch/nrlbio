#! /usr/bin/python
'''Sorts bed/gff file with a respect to strandness'''
import argparse
import sys;
import os;
from collections import defaultdict;

from pybedtools import BedTool;



parser = argparse.ArgumentParser(description='Sorts bed/gff file with a respect to strandness');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the bed/gff file to be sorted");
args = parser.parse_args();


sortedbed = BedTool(args.path).sort();
if(len(sortedbed) == 0):
	sys.stderr.write("input file is empty\n");
	sys.exit();
	

container = defaultdict(list)
container[sortedbed[0].strand].append(sortedbed[0])
chrom = sortedbed[0].chrom;

for interval in sortedbed[1:]:
	if(chrom == interval.chrom):
		container[interval.strand].append(interval)
		
	else:
		for strand, intervals in container.iteritems():
			for si in intervals:
				sys.stdout.write(str(si))
		chrom = interval.chrom;
		container = defaultdict(list)
		container[interval.strand].append(interval)
		
else:		
	for strand, intervals in container.iteritems():
		for si in intervals:
			sys.stdout.write(str(si))		
