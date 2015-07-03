#! /usr/lib/python
'''Generates decoy fasta file on basis of given one'''

import argparse
import os
import sys
import collections
from multiprocessing import Pool

from Bio import SeqIO

from nrlbio.generators import generator_seq
from nrlbio.HMM import MarkovChain

parser = argparse.ArgumentParser(description='Generates decoy fasta file on basis of given one');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file");
parser.add_argument('-o', '--order', nargs = '?', default = 1, type = int, help = "Markov model order. Correspondes to shuffling of \'order\'-nucleotides");
parser.add_argument('-ll', '--line_length', nargs = '?', default = 50, type = int, help = "Max length of line in output fasta file. Makes it more readable");
parser.add_argument('-t', '--threads', nargs = '?', default = 4, type = int, help = "Number of threads");
args = parser.parse_args();		


def gen_model(seqrecord):
	return MarkovChain.from_string(str(seqrecord.seq.upper()), args.order)
	
for seqrecord in SeqIO.parse(args.path, 'fasta'):
	mc = gen_model(seqrecord)
	print ">%s" % seqrecord.name
	for s in mc.generate_string(args.line_length):
		print s;	

#sys.exit()


#pool = Pool(processes = args.threads)
#res = pool.imap(gen_model, SeqIO.parse(args.path, 'fasta'))


#for mc in res:
	#for s in mc.generate_string(args.line_length):
		#print s;
	
				
