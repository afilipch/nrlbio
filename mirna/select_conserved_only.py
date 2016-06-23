#! /usr/lib/python
'''Selects only those miRNAs which have at least one[or more] conserved counterparts in other species'''

import argparse
import sys


from nrlbio.mirna import mirfasta2conservation


parser = argparse.ArgumentParser(description='Selects only those miRNAs which have at least one[or more] conserved counterparts in other species');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the consrevation fasta file of miRNAs");
parser.add_argument('--num', nargs = '?', default = 1, type = int, help = "miRNA is selected if it is conserved in [num] of species");
parser.add_argument('--species', nargs = '+', default = [], type = str, help = "miRNA is selected if it is conserved among [species]. If not set, this filter does not apply. [species] should follow mirBase name convention [hsa, mmu, cel]");
args = parser.parse_args();

species = set(args.species)

consdict, refmir = mirfasta2conservation(args.path)
#print refmir
#print
#print consdict

for mirid, d in consdict.items():
	if( (not species or species.issubset(set(d.keys()))) and len(d)>=args.num ):
		print ">%s\n%s" % (mirid, refmir[mirid])