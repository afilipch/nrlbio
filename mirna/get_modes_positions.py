#! /usr/lib/python
'''Finds seed-matches positions on miRNA targets'''

import argparse
import os
import sys

from nrlbio.mirna import fasta2mirnas
from nrlbio.generators import generator_doublebed

parser = argparse.ArgumentParser(description='Finds seed-matches positions on miRNA targets');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the mirna:target interactions, doublegff format");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs, fasta format");
args = parser.parse_args();






#get dictionary of mirna.Mirna objects
mirdict = fasta2mirnas(args.mir);


for i1, i2 in generator_doublebed(args.path):
	dpos = mirdict[i1.chrom].find_fixed_positions(i2.attrs['seq'])
	for type_, positions in dpos.items():
		print "%s\t%s\t%s\t%s\t%s" %  (type_, i2.chrom, i2.strand, i1.chrom,"\t".join([str(x+i2.start) for x in positions]));
	print
