#! /usr/lib/python
'''Searches for seed-matches types in given fasta sequences'''

import argparse
import os
import sys

from Bio import SeqIO 
from collections import defaultdict, Counter

from nrlbio.mirna import fasta2mirnas, assign_expression, mirnas2families, find_family, MODES_ORDER


parser = argparse.ArgumentParser(description='Searches for seed-matches types in given fasta sequences');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sequence to look for seed matches, fasta format");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs, fasta format");
args = parser.parse_args();






#get dictionary of mirna.Mirna objects
mirdict = fasta2mirnas(args.mir);


print "%s\t%s\t%s" % ('seq_id', 'mirna_id', "\t".join(MODES_ORDER));
for seqrecord in SeqIO.parse(args.path, "fasta"):
	tseq = str(seqrecord.seq.upper())
	for mirid, mirna in mirdict.items():
		modes = mirna.find_fixed_types(tseq);
		if(any(modes.values())):
			print "%s\t%s\t%s" % (seqrecord.name, mirid, "\t".join([str(modes[x]) for x in MODES_ORDER]));
		



	
	
	


	
