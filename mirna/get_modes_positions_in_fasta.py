#! /usr/lib/python
'''Searches for seed-matches types in given fasta sequences and outputs their positions'''

import argparse
import os
import sys

from Bio import SeqIO 
from collections import defaultdict, Counter

from nrlbio.mirna import fasta2mirnas, assign_expression, mirnas2families, find_family, MODES_ORDER
from nrlbio.sequencetools import shuffle_string

parser = argparse.ArgumentParser(description='Searches for seed-matches types in given fasta sequences');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sequence to look for seed matches, fasta format");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs, fasta format");
parser.add_argument('--shuffle', nargs = '?', default=False, const=True, type = bool, help = "if set target seqences will be shuffled");
args = parser.parse_args();






#get dictionary of mirna.Mirna objects
mirdict = fasta2mirnas(args.mir);


print "%s\t%s\t%s\t%s" % ('#seq_id', 'mirna_id', 'seed_type', "positions");
for seqrecord in SeqIO.parse(args.path, "fasta"):
	
	if(args.shuffle):
		tseq = shuffle_string(str(seqrecord.seq.upper()));
	else:
		tseq = str(seqrecord.seq.upper())
	
	for mirid, mirna in mirdict.items():
		dpos = mirna.find_fixed_positions(tseq)
		for type_, positions in dpos.items():
			if(positions):
				print "%s\t%s\t%s\t%s" %  (seqrecord.name, mirid, type_, ",".join([str(x) for x in positions]));
		



