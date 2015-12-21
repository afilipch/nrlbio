#! /usr/lib/python
'''Finds identity of provided intervals'''

import argparse
import os
import sys

from Bio import SeqIO, pairwise2
from collections import defaultdict, namedtuple



parser = argparse.ArgumentParser(description='Finds identity of provided intervals');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to sequences, fasta format");
parser.add_argument('--windows', nargs = '?', required=True, type = str, help = "path to the similar regions");
parser.add_argument('--refs', nargs = '?', required=True, type = str, help = "path to fasta reference from the first task");
args = parser.parse_args();

match_award = 2 
mm_penalty = -2.8
gap_penalty = -2.9;
sys.stderr.write("match_award: %1.2f\nmm_penalty%1.2f\ngap_penalty %1.2f\n" % (match_award, mm_penalty, gap_penalty));

def get_start_al(seq):
	for c, s in enumerate(seq):
		if(s!='-'):
			return c;


def get_end_al(seq):
	for c, s in enumerate(list(seq)[::-1]):
		if(s!='-'):
			return len(seq) - c;
			

def get_postitions(alignment, seq):
	aseq = alignment[1][alignment[3]:alignment[4]].replace('-', '')
	start = seq.find(aseq);
	end = start + len(aseq);
	return start, end


def get_identity(alignment):
	s1 = alignment[0][alignment[3]:alignment[4]]
	s2 = alignment[1][alignment[3]:alignment[4]]
	identity = 0;
	total = len(s1);
	for c1, c2 in zip(s1, s2):
		if(c1==c2):
			identity +=1;
	start = get_start_al(alignment[1])
	end = start+len(s1.replace("-", ""))

	return identity, total, len(s1.replace("-", "")), len(s2.replace("-", "")), start, end

	
	
		
def get_best_identity(seq1, seq2, ma, mp, gp):
	aset = set()
	result = [];
	identities = []
	for factor in range(10):
		ffactor = -factor/5.0
		for a in pairwise2.align.localms(seq2, seq1, ma, mp+ffactor, gp+ffactor, gp+ffactor):
			start1, end1 = get_postitions(a, seq1)
			#score, start2, end2 = a[2:]
			identity, total, l2, l1, start2, end2 = get_identity(a);
			fixed_total = len(seq1) - l1 + total 
			identities.append((identity/float(fixed_total), start2, end2, fixed_total))

	return max(identities, key= lambda x: x[0]);






seqs = {}
for path in args.path:
	for seqrecord in SeqIO.parse(path, "fasta"):
		seqs[seqrecord.name] = str(seqrecord.seq.upper());
		
		
regions = defaultdict(list);
with open(args.windows) as f:
	f.next()
	for l in f:
		a = l.strip().split("\t")
		id1, start1, end1, id2, start2, end2 = a[0], int(a[1]), int(a[2]), a[3], int(a[4]), int(a[5])
		regions[0].append((id1, start1, end1, seqs[id1][start1: end1]))
		regions[1].append((id2, start2, end2, seqs[id2][start2: end2]))
		
refs = {}
for c, seqrecord in enumerate(SeqIO.parse(args.refs, 'fasta')):
	refs[c] = (seqrecord.name, str(seqrecord.seq.upper()))
	
for count, rs in regions.items():
	refid, refseq = refs[count]
	for region in rs:
		identity, start2, end2, alength = get_best_identity(region[3], refseq, match_award, mm_penalty, gap_penalty)
		print "%s\t%d\t%d\t%s\t%d\t%d\t%1.2f\t%d" % (region[0], region[1], region[2], refid, start2, end2, identity*100, alength)
		
			 
			 


