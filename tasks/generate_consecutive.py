#! /usr/lib/python
'''Looks for pairs of windows with certain sequence identity'''

import argparse
import os
import sys

from Bio import SeqIO, pairwise2
from collections import defaultdict, namedtuple



parser = argparse.ArgumentParser(description='Looks for pairs of windows with certain sequence identity');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sequences, fasta format");
parser.add_argument('--minlength', nargs = '?', default=100, type = int, help = "min length of the windows");
parser.add_argument('--step', nargs = '?', default=20, type = int, help = "length of step between consecutive windows");
parser.add_argument('--identity', nargs = '?', default=0.65, type = float, help = "min sequence identity allowed");
parser.add_argument('--gaptype', nargs = '?', default='gap11', choices=['gap00', 'gap10', 'gap01', 'gap11'], type = str, help = "gap type");
args = parser.parse_args();

Piece = namedtuple('Piece', ['seq', 'start', 'end'])

def alfunc(seq1, seq2, ma, mp, gp, ftype):
	if(ftype=='gap11'):
		return pairwise2.align.localms(seq1, seq2, ma, mp, gp, gp);
	if(ftype=='gap10'):
		return pairwise2.align.localmd(seq1, seq2, ma, mp, gp, gp, -1000, -1000);
	if(ftype=='gap01'):
		return pairwise2.align.localmd(seq1, seq2, ma, mp, -1000, -1000, gp, gp);
	if(ftype=='gap00'):
		return pairwise2.align.localms(seq1, seq2, ma, mp, -1000, -1000);

def generate_piece(seq, length, step):
	ls = len(seq)
	#print ls;
	for start in range(0, ls-length-step, step):
		yield Piece(seq[start:start+length], start, start+length)
	else:	
		yield Piece(seq[start+step:], start+step, ls)
		
def get_postitions(alignment, seq):
	aseq = alignment[1][alignment[3]:alignment[4]].replace('-', '')
	start = seq.find(aseq);
	end = start + len(aseq);
	return start, end


def get_identity(alignment):
	s1 = alignment[0][alignment[3]:alignment[4]]
	s2 = alignment[1][alignment[3]:alignment[4]]
	identity = 0;
	total = float(len(s1));
	for c1, c2 in zip(s1, s2):
		if(c1==c2):
			identity +=1;
	
	return identity/total, total
	
	
		
def align(seq, piece, ma, mp, gp):
	#print "%d\t%d" % (piece.start, piece.end)
	aset = set()
	for a in alfunc(seq, piece.seq, ma, mp, gp, args.gaptype):
		start1, end1 = get_postitions(a, piece.seq)
		score, start2, end2 = a[2:]
		identity, alength = get_identity(a);
		if((start1, end1, start2, end2) not in aset):
			print "%d\t%d\t%d\t%d\t%1.2f\t%1.2f\t%d" % (start1+piece.start, end1+piece.start, start2, end2, score, identity, alength);
			aset.add((start1, end1, start2, end2))
	#print "*"*120
		
		

match_award = 2 
mm_penalty = -1*match_award*args.identity/(1-args.identity+0.01);
gap_penalty = mm_penalty - 0.05;
###gap_penalty = -1000;

sys.stderr.write("match_award: %1.2f\nmm_penalty%1.2f\ngap_penalty %1.2f\n" % (match_award, mm_penalty, gap_penalty))


seqs = []
for seqrecord in SeqIO.parse(args.path, 'fasta'):
	seqs.append(str(seqrecord.seq.upper()))
	
for pcount, piece in enumerate(generate_piece(seqs[0], args.minlength, args.step)):
	 align(seqs[1], piece, match_award, mm_penalty, gap_penalty)
	 sys.stderr.write('%d\n' % (pcount+1))
	
	
#print 


