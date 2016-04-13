#! /usr/bin/python
'''Finds seed matches for provided miRNAs. If expression of miRNAS is provided, than only seeds for top expressed miRNA families may be detected'''

import argparse
import os
import sys

from Bio import SeqIO 
from collections import defaultdict, Counter

from nrlbio.mirna import fasta2mirnas, assign_expression, mirnas2families, find_family
#from nrlbio.pyplot_extension import histogram

parser = argparse.ArgumentParser(description='Finds seed matches for provided miRNAs. If expression of miRNAS is provided, than only seeds for top expressed miRNA families may be detected');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sequence to look for seed matches, fasta format");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs, fasta format");
parser.add_argument('--expr', nargs = '?', type = str, help = "path to the expression file (mirid expression tsv format)");
parser.add_argument('--top', nargs = '?', default=10, type = int, help = "take only [--top] expressed miRNA families");
parser.add_argument('--start', nargs = '?', default=1, type = int, help = "seed start position on miRNA, 0-based inclusive");
parser.add_argument('--end', nargs = '?', default=7, type = int, help = "seed end position on miRNA, 0-based exclusive");
args = parser.parse_args();


def get_matches(seq, families):
	fam2matches = {};
	top_matches = [];
	for c, fam in enumerate(families):
		nmatches = len(fam.find_matches(seq, overlap=False));
		if(nmatches):
			fam2matches[fam.name] = nmatches;
			top_matches.append((c, 1));
	return fam2matches, top_matches


	



#get dictionary of mirna.Mirna objects
mirdict = fasta2mirnas(args.mir, args.start, args.end);

if(args.expr):
	#augment mirnas with expression info
	assign_expression(args.expr, mirdict, sep="\t");
	#filter out non-expressed mirnas
	mirdict = dict([ (x[0],x[1]) for x in mirdict.items() if x[1].expression])
	
#convolute mirnas into families
families = mirnas2families(mirdict.values())

if(args.expr and args.top>0):
	#sort families according to the expression
	families.sort(key=lambda x: x.expression, reverse = True)
	#get only the top expressed families
	families = families[:args.top];


matches_per_seq = []
top_fam_matches = defaultdict(int)
for seqrecord in SeqIO.parse(args.path, "fasta"):
	tseq = str(seqrecord.seq.upper())
	fam2matches, top_matches = get_matches(tseq, families);
	matches_per_seq.append(sum(fam2matches.values()));
	for k, v in top_matches:
		top_fam_matches[k] += v;
		


if(args.expr):
	cumsum = 0;
	for k in range(1, args.top):
		cumsum+=top_fam_matches[k];
		print "%d\t%d" % (k, cumsum);
	
	
	


	
