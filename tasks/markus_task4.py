#! /usr/lib/python
'''Finds identity of provided intervals'''

import argparse
import os
import sys

from Bio import SeqIO, pairwise2
from collections import defaultdict, namedtuple



parser = argparse.ArgumentParser(description='Finds identity of provided intervals');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to fasta reference from the first task");
parser.add_argument('--windows', nargs = '?', required=True, type = str, help = "output of task 3");
parser.add_argument('--oal', nargs = '?', default=False, const=True, type = bool, help = "if set, alignment is output");
args = parser.parse_args();

match_award = 2 
mm_penalty = -2.3
gap_penalty = -2.4;
sys.stderr.write("match_award: %1.2f\nmm_penalty%1.2f\ngap_penalty %1.2f\n" % (match_award, mm_penalty, gap_penalty));

def get_start_al(seq):
	for c, s in enumerate(seq):
		if(s!='-'):
			return c;


def fix_test_alignment(alstr):
	m = alstr.split("\n")
	s1 = [];
	s2 = [];
	#slines = [];
	for c1,cl, c2 in zip(m[0], m[1], m[2]):
		#if(cl=="|"):
			#slines.append(cl);
			
		if(c1 == c2):
			s1.append(c1);
			s2.append(c2);
		else:
			s1.append(c1.lower());
			s2.append(c2.lower());
		
	m[0] = "".join(s1);
	#m[1] = "".join(slines);
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
	for factor in range(40):
		ffactor = -factor/10.0
		for a in pairwise2.align.localms(seq2, seq1, ma, mp+ffactor, gp+ffactor, gp+ffactor):
			start1, end1 = get_postitions(a, seq1)
			identity, total, l2, l1, start2, end2 = get_identity(a);
			fixed_total = len(seq1) - l1 + total 
			identities.append((identity/float(fixed_total), start2, end2, fixed_total, fix_test_alignment(pairwise2.format_alignment(*a))))

	return max(identities, key= lambda x: x[0]);






seqs = {}
for path in args.path:
	for seqrecord in SeqIO.parse(path, "fasta"):
		seqs[seqrecord.name] = str(seqrecord.seq.upper());
		
		
regions = []
with open(args.windows) as f:
	for l in f:
		a = l.strip().split("\t")
		id1, start1, end1, id2, start2, end2 = a[0], int(a[1]), int(a[2]), a[3], int(a[4]), int(a[5])
		regions.append((id2, start2, end2, seqs[id2][start2: end2]))

hpv16 =  regions[:len(regions)/2]
hpv18 =  regions[len(regions)/2:]

for r16, r18 in zip(hpv16, hpv18):
	identity, start2, end2, alength, astr = get_best_identity(r16[3], r18[3], match_award, mm_penalty, gap_penalty)
	print "%s\t%d\t%d\t%s\t%d\t%d\t%1.2f\t%d" % (r16[0], r16[1], r16[2], r18[0], r18[1], r18[2], identity*100, alength)
	if(args.oal):
		print astr
		print
		
			 
			 


