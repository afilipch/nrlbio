#! /usr/lib/python
'''Finds identity of provided intervals'''

import argparse
import os
import sys

from Bio import SeqIO, pairwise2
from collections import defaultdict, namedtuple



parser = argparse.ArgumentParser(description='Finds identity of provided intervals');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to sequences, fasta format");
parser.add_argument('--windows', nargs = 2, required=True, type = str, help = "path to the similar regions");
args = parser.parse_args();

match_award = 2 
mm_penalty = -2.5
gap_penalty = -2.6;
sys.stderr.write("match_award: %1.2f\nmm_penalty%1.2f\ngap_penalty %1.2f\n" % (match_award, mm_penalty, gap_penalty));



def get_identity(alignment):
	s1 = alignment[0][alignment[3]:alignment[4]]
	s2 = alignment[1][alignment[3]:alignment[4]]
	identity = 0;
	total = float(len(s1));
	for c1, c2 in zip(s1, s2):
		if(c1==c2):
			identity +=1;
	
	return identity, total, len(s1.replace("-", "")), len(s2.replace("-", ""))
	
	
		
def get_best_identity(seq1, seq2, ma, mp, gp):
	aset = set()
	result = [];
	identities = [];
	for n in range(0,50):
		add = -n/10.0
		for a in pairwise2.align.localms(seq1, seq2, ma, mp+add, gp+add, gp+add):
			identity, total, l1, l2 = get_identity(a);
			fixed_total = min(len(seq1)-l1+total, len(seq2)-l2+total)
			identities.append(identity/fixed_total);
	return max(identities);

#regions1 = [];
#with open(args.windows[0]) as f:
	#f.next()
	#for l in f:
		#a = l.strip().split("\t")
		#regions1.append((a[0], int(a[1]), int(a[2])))
		#regions1.append((a[3], int(a[4]), int(a[5])))
	
	
#regions2 = [];	
#with open(args.windows[1]) as f:
	#f.next()
	#for l in f:
		#a = l.strip().split("\t")
		#regions2.append((a[0], int(a[1]), int(a[2])))
		#regions2.append((a[3], int(a[4]), int(a[5])))		
		

seqs = {}
for path in args.path:
	for seqrecord in SeqIO.parse(path, "fasta"):
		seqs[seqrecord.name] = str(seqrecord.seq.upper());
		
		
regions1 = [];
with open(args.windows[0]) as f:
	f.next()
	for l in f:
		a = l.strip().split("\t")
		id1, start1, end1, id2, start2, end2 = a[0], int(a[1]), int(a[2]), a[3], int(a[4]), int(a[5])
		regions1.append((id1, start1, end1, seqs[id1][start1: end1]))
		regions1.append((id2, start2, end2, seqs[id2][start2: end2]))
		
		
regions2 = [];
with open(args.windows[1]) as f:
	f.next()
	for l in f:
		a = l.strip().split("\t")
		id1, start1, end1, id2, start2, end2 = a[0], int(a[1]), int(a[2]), a[3], int(a[4]), int(a[5])
		regions2.append((id1, start1, end1, seqs[id1][start1: end1]))
		regions2.append((id2, start2, end2, seqs[id2][start2: end2]))

print "\t".join(("ID_1", "Start_1", "End_1", "ID_2", "Start_2", "End_2", "Identity"))		
for r1 in regions1:
	for r2 in regions2:
		identity =get_best_identity(r1[3], r2[3], match_award, mm_penalty, gap_penalty)
		print "%s\t%d\t%d\t%s\t%d\t%d\t%1.2f" % (r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], identity)
		













	
