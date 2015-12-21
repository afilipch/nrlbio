#! /usr/lib/python
'''Checks if the seed match covers the splice junction'''

import argparse
import os
import sys
from collections import defaultdict, Counter
import random

from Bio import SeqIO

from nrlbio.mirna import fasta2mirnas
from nrlbio.generators import generator_doublebed


parser = argparse.ArgumentParser(description='Assignes binding mode for mirna:target interaction');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the mirna:target interactions, doublegff format");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs, fasta format");
parser.add_argument('--junctions', nargs = '?', required=True, type = str, help = "path to the junctions sequences, fasta format");
parser.add_argument('--length', nargs = '?', required = True, type = int, help = "length of splice junction")
parser.add_argument('--candidates', nargs = '?', required = True, type = str, help = "Path to the output file for interactions with seed match covering splice junction")
args = parser.parse_args();

bottom = args.length/2 - 7 
upper = args.length/2+1


def evaluate_pos(pos):
	if(pos>-1):
		return range(pos,pos+6)
	else:
		return None;
	

		
		
		
mirnas = fasta2mirnas(args.mir)

junction2mirid = {}
for i1, i2 in generator_doublebed(args.path):
	junction2mirid[i2.chrom] = i1.chrom


coverage = [];
candidates = [];

for seqrecord in SeqIO.parse(args.junctions, 'fasta'):
	mirid = junction2mirid.get(seqrecord.name, None)
	if(mirid):
		mirna = mirnas[junction2mirid[seqrecord.name]]
		#mirna = random.choice(mirnas.values())
		seq = str(seqrecord.seq.upper())
		pos = seq.find(mirna.match)
		positions = evaluate_pos(pos)
		if(positions):
			coverage.extend(positions);
		if(pos > bottom and pos < upper):
			candidates.append((mirid, seqrecord.name));
			
			
	
coverage = Counter(coverage);
	
for i in range(args.length):
	print "%d\t%d" % (i, coverage[i])
	
	
with open(args.candidates, 'w') as f:
	for i1, i2 in generator_doublebed(args.path):
		if((i1.chrom, i2.chrom) in candidates):
			f.write(str(i1))
			f.write(str(i2))
