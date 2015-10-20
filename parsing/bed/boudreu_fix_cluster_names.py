#! /usr/bin/python
'''Converts boudreu bed intervals names into ones compatible with assign_coordinates.py script'''
import argparse
import sys;

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts boudreu bed intervals names into ones compatible with assign_coordinates.py script');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the bed/gff file");
args = parser.parse_args();

for i in BedTool(args.path):
	a = i.chrom.split("|")
	if(len(a)==2):
		i.chrom = "|".join(a[1].split(":"));
	sys.stdout.write(str(i))