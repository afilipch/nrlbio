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
args = parser.parse_args();		

nucleotides = set('ACTG')

l = [];	
for seqrecord in SeqIO.parse(args.path, 'fasta'):
	print ">random_%s" % seqrecord.name
	
	for n in seqrecord.seq.upper():
		if(n=="N"):
			l.append(n);
		elif(random.random()>args.rate):
			l.append(n);
		else:
			m = random.choice(list(nucleotides-set(n)));
			l.append(m);
			
		if(len(l) == args.line_length):
			print "".join(l)
			l[:] = []

			
	else:
		if(l):
			print "".join(l);
		