#! /usr/bin/python
'''Selects interactions only for mirna ids provided''' 
import argparse
import sys;

from nrlbio.generators import generator_doublebed



parser = argparse.ArgumentParser(description='Selects interactions only for mirna ids provided');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the mirna:target interactions to select from, double gff");
parser.add_argument('--mirids', nargs = '+', type = str, help = "miRNA ids to select");
parser.add_argument('--inverse', nargs = '?', default=False, const=True, type = bool, help = "inverse selection. Select all interactions with miRNAs not for [--mirids]");
args = parser.parse_args();

mirids = set(args.mirids)

if(args.inverse):
	for i1, i2 in generator_doublebed(args.path):
		if(i1.chrom not in mirids):
			sys.stdout.write(str(i1))
			sys.stdout.write(str(i2))
else:		
	for i1, i2 in generator_doublebed(args.path):
		if(i1.chrom in mirids):
			sys.stdout.write(str(i1))
			sys.stdout.write(str(i2))
		