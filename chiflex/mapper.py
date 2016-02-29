#! /usr/lib/python
'''Custom mapper useful to map to miRNAs or other small reference'''

import argparse
import os
import sys
import timeit
from collections import defaultdict, namedtuple


from Bio import SeqIO, pairwise2



parser = argparse.ArgumentParser(description='Finds identity of provided intervals');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sequences, fasta/fastq format");
parser.add_argument('--reference', nargs = '?', required=True, type = str, help = "path to the reference in fasta format");
parser.add_argument('--rtype', nargs = '?', default='fastq', choices=['fasta', 'fastq'], type = str, help = "Type of input reads, fasta or fastq");
parser.add_argument('--plength', nargs = '?', default=13, type = int, help = "Only first [plength] nucleotides will be used for anchoring");
parser.add_argument('--minscore', nargs = '?', default=26, type = int, help = "Min aligned score of a hit to be reprted");
args = parser.parse_args();


def find_matches(seq, reference, anchors, plength, minscore):
	fseq = seq[:plength*2]
	bestmatches = [];	
	for name in anchors[seq[:plength]]:
		matches = pairwise2.align.localms(fseq, reference[name], 2, -5, -4, -6)
		if(matches and matches[0][2]>minscore):
			bestmatches.append((name, matches[0]));
	return bestmatches;




reference = {}
for seqrecord in SeqIO.parse(args.reference, 'fasta'):
	reference[seqrecord.name] = str(seqrecord.seq.upper());
	
anchors = defaultdict(list);
for name, seq in reference.items():
	anchors[seq[:args.plength]].append(name);
	
	

start_time = timeit.default_timer()
for count, seqrecord in enumerate(SeqIO.parse(args.path, args.rtype)):
	seq = str(seqrecord.seq);
	matches = find_matches(seq, reference, anchors, args.plength, args.minscore);
	print len(matches);
	if((count+1) % 1000 == 0):
		sys.stderr.write("%d reads mapped in %1.1f seconds\n" % (count+1, timeit.default_timer() - start_time))
		
		
		
		