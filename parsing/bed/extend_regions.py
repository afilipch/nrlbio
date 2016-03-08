#! /usr/bin/python
'''Extends all regions in the provided bed/gff file to a fixed value''' 
import argparse
import sys;

from pybedtools import BedTool



parser = argparse.ArgumentParser(description='Extends target region of  miRNA:target interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mirna:target interactions to select from, double gff/bed");
parser.add_argument('-l', '--left', nargs = '?', required = True, type = int, help = "Left extension in nucleotides");
parser.add_argument('-r', '--right', nargs = '?', required = True, type = int, help = "Right extension in nucleotides");
parser.add_argument('-s', '--strand', nargs = '?', default = False, const = True, type = bool, help = "If set left and right extensions are swapped for the minus strand");
args = parser.parse_args();

	
for interval in BedTool(args.path):
	if(args.strand and interval.strand == '-'):
		left = args.right;
		right = args.left;
	else:
		right = args.right;
		left = args.left;
		
	interval.start = max(interval.start-left, 0)
	interval.end = interval.end+right;
	
	sys.stdout.write(str(interval))