#! /usr/lib/python
'''Generates decoy fasta file on basis of given one'''

import argparse
import os;
import sys

from Bio import SeqIO

from nrlbio.sequencetools import split2chunks;
import nrlbio.HMM as hmm;

parser = argparse.ArgumentParser(description='Generates decoy fasta file on basis of given one');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file");
parser.add_argument('-m', '--multimodel', nargs = '?', default = False, const = True, type = bool, help = "Use set of Markov models to generate decoy sequence");
parser.add_argument('-o', '--order', nargs = '?', default = 1, type = int, help = "Markov model order. Correspondes to shuffling of \'order\'-nucleotides");
parser.add_argument('-w', '--window', nargs = '?', default = 400, type = int, help = "Used in multimodel regime. Each fasta sequence will be split into sequences of length \'window\', for each of them separate Markov Model will be generated");
parser.add_argument('-ll', '--line_length', nargs = '?', default = 100, type = int, help = "max length of line in output fasta file. Makes it more readable");
args = parser.parse_args();

#sys.stderr.write("%s\n" % args.multimodel)
	

if(args.multimodel):
	for seq_record in SeqIO.parse(args.path, "fasta"):
		print ">" + "_".join(["random", seq_record.id]);
		
		mm = hmm.Meta_Markov_Model.from_sequence(str(seq_record.seq.upper()), args.order, args.window)
		
		seq = '' 
		for s in mm.generate_string(len(seq_record.seq)):
			seq += s;
			nchunks = len(seq)/args.line_length
			if(nchunks):
				for chunk in split2chunks(seq[:nchunks*args.window], args.line_length):
					print chunk;
				seq = seq[nchunks*args.window:]
			else:
				pass;
		if(seq):		
			print seq;
				
else:
	for seq_record in SeqIO.parse(args.path, "fasta"):
		print "_".join([">" + seq_record.id, "random"]);
		
		mm = hmm.Markov_Model.from_sequence(str(seq_record.seq.upper()), args.order)
		seq = mm.generate_string(len(seq_record.seq));
		
		for chunk in split2chunks(seq, args.line_length):
			print chunk;	
	