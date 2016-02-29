#! /usr/bin/python
'''Outputs targets of miRNA''' 
import argparse
import sys;

from nrlbio.generators import generator_doublebed



parser = argparse.ArgumentParser(description='Outputs targets of miRNA');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the mirna:target interactions to select from, double gff/bed");
args = parser.parse_args();

	
for i1, i2 in generator_doublebed(args.path):
	i2.name = i1.chrom
	sys.stdout.write(str(i2))