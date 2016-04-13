#! /usr/bin/python
'''Checks miRNA star/mature choice determinants''' 
import argparse
import sys;
import math;
from collections import defaultdict

from Bio import SeqIO;
import jinja2

from nrlbio.rnahybrid import get_rnahybrid


parser = argparse.ArgumentParser(description='Checks miRNA star/mature choice determinantsg');
parser.add_argument('--hairpins', nargs = '?', required = True, type = str, help = "Path to the miRNA hairpins, fasta file")
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the miRNAs, fasta format");
parser.add_argument('--expression', nargs = '?', required = True, type = str, help = "Path to the miRNA expression, tsv file");
args = parser.parse_args();


def get2score(pattern, energy):
	return (sum(pattern[1:8]) - sum(pattern[9:12]) + sum(pattern[12:]) + math.log(-1*energy, 2))/float(len(pattern))


hairpins = {}
for seqrecord in SeqIO.parse(args.hairpins, 'fasta'):
	name = '-'.join(seqrecord.id.split('-')[:3]).replace('mir', 'miR')
	seq = str(seqrecord.seq.upper()).replace('U', 'T')
	l = len(seq)/2
	hairpins[name] = (seq[:l], seq[l:])
	
mirna5 = {}
mirna3 = {}

for seqrecord in SeqIO.parse(args.mir, 'fasta'):
	a = seqrecord.id.split('-')
	name = '-'.join(a[:3])
	seq = str(seqrecord.seq.upper()).replace('U', 'T')
	if(a[-1]=='3p'):
		mirna3[name] =  seq
	elif(a[-1]=='5p'):
		mirna5[name] =  seq
	else:
		pass;
	
	
#print len(set(mirna5.keys()) & set(hairpins.keys()) )

	
dscores = {}

for mirid, (seq5, seq3) in hairpins.items():
	mir5 = mirna5.get(mirid)
	mir3 = mirna3.get(mirid)
	
	if(mir5 and mir3):
		energy3, pattern3 = get_rnahybrid(seq5, mir3);
		score3 = get2score(pattern3, energy3);
		
		energy5, pattern5 = get_rnahybrid(seq3, mir5);
		score5 = get2score(pattern5, energy5);
		
		dscores[mirid] = (score5, score3)
	
	
	
expression5 = defaultdict(float)
expression3 = defaultdict(float)

with open(args.expression) as f:
	f.next();
	for l in f:
		arr = l.strip().split("\t")
		a = arr[0].split("-")
		name = '-'.join(a[:3])
		if(a[-1]=='3p'):
			expression3[name] +=  sum([float(x) for x in arr[1:]])
		elif(a[-1]=='5p'):
			expression5[name] +=  sum([float(x) for x in arr[1:]])
		else:
			pass;
		
expression = {}
for k, e5 in expression5.items():
	e3 = expression3[k]
	if(e5 and e3):
		expression[k] = (e5, e3)
	
	
for mirid, (ds5, ds3) in dscores.items():
	e5, e3 = expression.get(mirid, 'nn')
	if(e5!='n'):
		print "%1.3f\t%d" % ((ds5-ds3), (e5-e3))
		
		

