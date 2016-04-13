#! /usr/bin/python
'''Flatten each interaction with ambiguous mirna into separate interactions with each of the ambiguous mirna'''

import argparse
import os
import sys

from Bio import SeqIO 
from collections import defaultdict

from nrlbio.generators import generator_doublebed


parser = argparse.ArgumentParser(description='Flatten each interaction with ambiguous mirna into separate interactions with each of the ambiguous mirna');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the miRNA:target interactions");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs to assign sequences to new interactions. Fasta format");
args = parser.parse_args();

id2seq = {}

for seqrecord in SeqIO.parse(args.mir, 'fasta'):
	id2seq[seqrecord.name] = str(seqrecord.seq.upper()).replace('U', 'T')
	
	

for i1, i2 in generator_doublebed(args.path):
	a = i1.chrom.split(",")
	if(len(a) > 1):
		for mirid in a:
			i1.attrs['seq'] = id2seq[mirid]
			sys.stdout.write('%s%s' % (i1, i2))
	else:
		sys.stdout.write('%s%s' % (i1, i2))
