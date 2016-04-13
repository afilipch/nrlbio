#! /usr/bin/python
'''Checks provided sequences for extendet complemetarity''' 
import argparse
import sys;
from itertools import product

from Bio import pairwise2, SeqIO
from Bio.Seq import reverse_complement



parser = argparse.ArgumentParser(description='Checks provided sequences for extendet complemetarity');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "Path to the sequences. All sequences from the first file will be aligned against all from the second one. Two files in fasta format");
args = parser.parse_args();

seqs1 = []
for seqrecord in SeqIO.parse(args.path[0], 'fasta'):
	seqs1.append((seqrecord.id, str(seqrecord.seq.upper()).replace('U', 'T')))

seqs2 = []
for seqrecord in SeqIO.parse(args.path[1], 'fasta'):
	seqs2.append((seqrecord.id, reverse_complement(str(seqrecord.seq.upper()).replace('U', 'T'))))
	
for (id1, seq1), (id2, seq2) in product(seqs1, seqs2):
	a = pairwise2.align.localms(seq1, seq2, 2, -5, -4, -3)
	print
	print id1, id2
	print a[0]
	print 
	print "_"*140
	







