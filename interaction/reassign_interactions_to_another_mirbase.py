#! /usr/bin/python
'''reassign miRNA ids to a new ones from a given mirbase'''

import argparse
import os
import sys

from Bio import SeqIO 
from collections import defaultdict

from nrlbio.generators import generator_doublebed


parser = argparse.ArgumentParser(description='reassign miRNA ids to a new ones from a given mirbase');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the miRNA:target interactions");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs. New Id's will be taken from this file. Fasta format");
args = parser.parse_args();

seq2id = {}

for seqrecord in SeqIO.parse(args.mir, 'fasta'):
	seq2id[str(seqrecord.seq.upper()).replace('U', 'T')] = seqrecord.name
	
#for k, v in seq2id.items():
	#if(len(v)>1):
		#print k, v
		#print
	

skipped = 0
for i1, i2 in generator_doublebed(args.path):
	a = i1.chrom
	new= seq2id.get(i1.attrs['seq'], None)
	if(new):
		i1.chrom = new
		sys.stdout.write('%s%s' % (i1, i2))
	else:
		for seq, mirid in seq2id.items():
			if( seq.startswith(i1.attrs['seq']) or i1.attrs['seq'].startswith(seq) ):
				i1.chrom = mirid
				i1.attrs['seq'] = seq
				sys.stdout.write('%s%s' % (i1, i2))
				break
		else:
			sys.stderr.write('%s\t%s\n' % (i1.chrom, i1.attrs['seq']))
			skipped += 1;
		
sys.stderr.write('%d interactions were skipped because of missing counterparts\n' % skipped)