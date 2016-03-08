#! /usr/lib/python
'''Checks if UTR's of differential length carry extra miRNA targets'''

import argparse
import os
import sys
from collections import defaultdict



parser = argparse.ArgumentParser(description='Checks if UTR\'s of differential length carry extra miRNA targets');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to miRNA's binding sites counts in UTR");
parser.add_argument('--merged', nargs = '?', default=False, const=True, type = bool, help = "If set, then binding sites will be merged for different miRNA");
##parser.add_argument('-sc', '--set_control', nargs = '?', default = False, const=True, type = int, help = "If set, background(shuffled target sequence) seed match mode will be assigned");
args = parser.parse_args();

llimit = 20


def seed_score(counts):
	return counts[0]*3 + counts[1]*2 + counts[2]*2 + counts[3] + counts[4]*2 + counts[5]*3 

def broad_score(counts):
	m29a, m28a, m27a, m29, m28, m27, m38, mm28 = counts
	return 0.1*(m29a + m29 + m28a) + 0.2*(m28+m27a) + m27 + 0.25*m38 + 0.25*mm28

def compare_mirna(utrs, utr_name, scorefun=seed_score):
	iso1, iso2 = utrs.values();
	l1, l2 = float(iso1[0][0]), float(iso2[0][0])
	if(l1 < iso2[0][0]):
		iso1, iso2 = iso2, iso1
		l1, l2 = l2, l1
		
	iso1 = dict([ (x[1], scorefun(x[2])) for x in iso1 ])
	iso2 = dict([ (x[1], scorefun(x[2])) for x in iso2 ])
	
	for mirid, score1 in iso1.items():
		score2 = iso2.get(mirid, 0)
		if(score1 or score2):
			if(llimit<l1-l2):
				print "%s\t%s\t%1.4f\t%1.4f\t%1.4f" % (utr_name, mirid, score1/l1, score2/l2, (score1-score2)/(l1-l2))
			
			
def compare_merged(utrs, utr_name, scorefun=broad_score):
	iso1, iso2 = utrs.values();
	l1, l2 = float(iso1[0][0]), float(iso2[0][0])
	if(l1 < iso2[0][0]):
		iso1, iso2 = iso2, iso1
		l1, l2 = l2, l1
		
	score1 = sum([scorefun(x[2]) for x in iso1[1:]])
	score2 = sum([scorefun(x[2]) for x in iso2[1:] ])

	if(score1 or score2):
		if(llimit < (l1-l2) and (score1-score2)):
			s1 = score1/l1
			s2 = score2/l2
			sgain = (score1-score2)/(l1-l2)
			print "%s\t%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % (utr_name, 'merged', s1, s2, sgain, sgain-s2)
			

if(args.merged):
	compare = compare_merged
else:
	compare = compare_mirna

meta_utrs = defaultdict(lambda: defaultdict(list))

with open(args.path) as f:
	f.next()
	for l in f:
		a = l.strip().split("\t");
		utr_name, utr_isoform = a[0].split("|");
		meta_utrs[utr_name][utr_isoform].append(( int(a[1]), a[2], tuple([int(x) for x in a[3:]]) ));
		
		

for utr_name, utrs in meta_utrs.items():
	compare(utrs, utr_name);