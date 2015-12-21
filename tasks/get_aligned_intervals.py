#! /usr/lib/python
'''Finds exact alignments on basis of alignment of pieces'''

import argparse
import os
import sys

from Bio import SeqIO, pairwise2
from collections import defaultdict, namedtuple



parser = argparse.ArgumentParser(description='Finds exact alignments on basis of alignment of pieces');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sequences, fasta format");
parser.add_argument('--pieces', nargs = '?', required=True, type = str, help = "path to aligned pieces");

parser.add_argument('--minlength', nargs = '?', default=40, type = int, help = "min length of the windows");
parser.add_argument('--identity', nargs = '?', default=0.75, type = float, help = "min sequence identity allowed");
parser.add_argument('--gaptype', nargs = '?', default='gap11', choices=['gap00', 'gap10', 'gap01', 'gap11'], type = str, help = "gap type");
parser.add_argument('--oal', nargs = '?', default=False, const=True, type = bool, help = "if set, alignment is output");
args = parser.parse_args();

match_award = 2 
mm_penalty = -1*match_award*args.identity/(1-args.identity+0.01);
gap_penalty = mm_penalty - 0.05;
sys.stderr.write("match_award: %1.2f\nmm_penalty%1.2f\ngap_penalty %1.2f\n" % (match_award, mm_penalty, gap_penalty))

def alfunc(seq1, seq2, ma, mp, gp, ftype):
	if(ftype=='gap11'):
		return pairwise2.align.localms(seq1, seq2, ma, mp, gp, gp);
	if(ftype=='gap10'):
		return pairwise2.align.localmd(seq1, seq2, ma, mp, gp, gp, -1000, -1000);
	if(ftype=='gap01'):
		return pairwise2.align.localmd(seq1, seq2, ma, mp, -1000, -1000, gp, gp);
	if(ftype=='gap00'):
		return pairwise2.align.localms(seq1, seq2, ma, mp, -1000, -1000);


def fix_test_alignment(alstr):
	m = alstr.split("\n")
	s1 = [];
	s2 = [];
	slines = [];
	for c1,cl, c2 in zip(m[0], m[1], m[2]):
		if(cl=="|"):
			slines.append(cl);
			
			if(c1 == c2):
				s1.append(c1);
				s2.append(c2);
			else:
				s1.append(c1.lower());
				s2.append(c2.lower());
			
	m[0] = "".join(s1);
	m[1] = "".join(slines);
	m[2] = "".join(s2);
	
	return "\n".join(m[:4])

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
	
	return identity/total, total, len(s1.replace("-", "")), len(s2.replace("-", ""))
	
	
		
def align(seq1, seq2, window, ma, mp, gp, id1, id2, oal):
	aset = set()
	result = [];
	for a in alfunc(seq2, seq1, ma, mp, gp, args.gaptype):
		start1, end1 = get_postitions(a, seq1)
		score, start2, end2 = a[2:]
		identity, alength, l1, l2 = get_identity(a);
		if(identity>=0.75 and l1>=40 and l2>=40):
			if((start1, end1, start2, end2) not in aset):
				print "%s\t%d\t%d\t%s\t%d\t%d\t%1.2f\t%d" % (id1, start1+window.start1, start1+window.start1+l1, id2, start2+window.start2, start2+window.start2+l2, identity*100, a[4] - a[3]);
				if(oal):
					print fix_test_alignment(pairwise2.format_alignment(*a))
					print 
				#print a[0]
				#print a[1]
				aset.add((start1, end1, start2, end2))







AlignedPiece = namedtuple('Piece', ['start1', 'end1', 'start2', 'end2'])

def check_window(piece, window):
	d1 = min(piece.end1, window.end1) - max(piece.start1, window.start1)
	d2 = min(piece.end2, window.end2) - max(piece.start2, window.start2)
	return d1>=0 and d2>=0;

pieces = [];
with open(args.pieces) as f:
	for l in f:
		a = [int(x) for x in l.strip().split("\t")[:4]]
		pieces.append(AlignedPiece(*a))
		
pieces.sort(key =  lambda x: x.start1)

windows = []

for piece in pieces:
	for index, window in enumerate(windows):
		if(check_window(piece, window)):
			windows[index] = AlignedPiece(min(piece.start1, window.start1), max(piece.end1, window.end1), min(piece.start2, window.start2), max(piece.end2, window.end2))
			break;
	else:
		windows.append(piece)
			
			




seqs = []
ids = []
for seqrecord in SeqIO.parse(args.path, 'fasta'):
	seqs.append(str(seqrecord.seq.upper()))
	ids.append(seqrecord.name);
	
print "\t".join(("ID_1", "Start_1", "End_1", "ID_2", "Start_2", "End_2", "Identity", "length_alignment"))
	
for window in windows:
	seq1 = seqs[0][window.start1:window.end1]
	seq2 = seqs[1][window.start2:window.end2]
	align(seq1, seq2, window, match_award, mm_penalty, gap_penalty, ids[0], ids[1], args.oal)
	
	
#windows.sort(key =  lambda x: x.start1)

#for window in windows:
	#print window


#for piece in pieces:
	#print piece