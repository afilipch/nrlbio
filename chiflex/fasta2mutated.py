#! /usr/lib/python
'''Generates fasta file with mutations on basis of given one'''

import argparse
import os
import sys
import random


from Bio import SeqIO



parser = argparse.ArgumentParser(description='Generates fasta file with mutations on basis of given one');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file");
parser.add_argument('-r', '--rate', nargs = '?', default = 0.1, type = float, help = "rate of mutations, float [0,1]");
parser.add_argument('-ll', '--line_length', nargs = '?', default = 50, type = int, help = "Max length of line in output fasta file. Makes it more readable");
parser.add_argument('--local', nargs = '?', const=True, default = False, type = float, help = "If set, mutations will be distributed uniformaly, with slightly increasing rate(edge effect). This option should be activated in case of short sequences(miRNAs)");
args = parser.parse_args();

if(args.rate>1 or args.rate<0):
	sys.exit('Error:: Mutation rate has to be in a range [0,1]. It was set to %1.2f' % args.rate)
else:
	rate = args.rate

each = int(1/rate)


nucleotides = set('ACTG')
mdict = dict([ (x, tuple(nucleotides - set(x))) for x in 'ACTG']);


def mutate_random(seqrecord, each):
	l = [];
	for n in seqrecord.seq.upper():
		if(n=="N"):
			l.append(n);
		else:	
			if(random.random()>rate):
				l.append(n);
			else:
				m = random.choice(mdict[n]);
				l.append(m);
			
		if(len(l) == args.line_length):
			print "".join(l)
			l[:] = []
			
	if(l):
		print "".join(l);


def mutate_local(seqrecord, each):
	piece = []
	l = [];
	for n in seqrecord.seq.upper():
		piece.append(n);
		if(len(piece)==each):
			rindex = random.randrange(0,len(piece))
			piece[rindex] = random.choice(mdict[piece[rindex]]);
			l.extend(piece);
			piece[:] = []
		
	if(piece):
		rindex = random.randrange(0,len(piece))
		piece[rindex] = random.choice(mdict[piece[rindex]]);
		l.extend(piece);
		
	print "".join(l);
			
if(args.local):
	mutate = mutate_local;
else:
	mutate = mutate_random

for seqrecord in SeqIO.parse(args.path, 'fasta'):
	print ">random_%s" % seqrecord.name
	mutate(seqrecord, each)



	