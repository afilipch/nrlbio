#! /usr/bin/python
'''Extends target region of  miRNA:target interactions''' 
import argparse
import sys;

from nrlbio.generators import generator_doublebed



parser = argparse.ArgumentParser(description='Extends target region of  miRNA:target interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mirna:target interactions to select from, double gff/bed");
parser.add_argument('-l', '--left', nargs = '?', required = True, type = int, help = "Left extension in nucleotides");
parser.add_argument('-r', '--right', nargs = '?', required = True, type = int, help = "Right extension in nucleotides");
args = parser.parse_args();

	
for i1, i2 in generator_doublebed(args.path):
	i2.start = max(i2.start-args.left, 0)
	i2.end = i2.end+args.right;
	sys.stdout.write(str(i1))
	sys.stdout.write(str(i2))