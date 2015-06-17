#! /usr/lib/python
'''Generates decoy fasta file on basis of given one'''

import argparse
import os
import sys
import collections
from multiprocessing import Pool

from Bio import SeqIO

from nrlbio.sequencetools import split2chunks;
from nrlbio.HMM import MarkovChain, MultiMarkov, SlidingMarkovChain, GrowingMarkovChain, markov_difference, find_switch

parser = argparse.ArgumentParser(description='Generates decoy fasta file on basis of given one');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file");
parser.add_argument('-m', '--multimodel', nargs = '?', default = False, const = True, type = bool, help = "Use set of Markov models to generate decoy sequence");
parser.add_argument('-o', '--order', nargs = '?', default = 1, type = int, help = "Markov model order. Correspondes to shuffling of \'order\'-nucleotides");
parser.add_argument('-w', '--window', nargs = '?', default = 400, type = int, help = "Can be used in multimodel regime. Each fasta sequence will be split into sequences of length \'window\', for each of them separate Markov Model will be generated");
parser.add_argument('-ll', '--line_length', nargs = '?', default = 100, type = int, help = "Max length of line in output fasta file. Makes it more readable");
args = parser.parse_args();


		

#models = [];
#step = 50000
for seqrecord in SeqIO.parse(args.path, "fasta"):
	#mm = MultiMarkov.from_seqrecord(seqrecord, args.order)
	#mm.serialize("%s.yml" % seqrecord.name)
	#mm = MultiMarkov.deserialize("%s.yml" % seqrecord.name)
	#for i, model in enumerate(mm.models):
		#print "model length: %d" % model.length
		#for fr, d in model.transitions.items():
			#print "\nfrom %s\ttotal: %d" % (fr, sum(d.values()))
			#for t, p in d.items():
				#print "\tto: %s\tprob: %d" % (t, p)
		#print "\n\n___________________________________________________________________________________________________________________"
sys.exit();	
						

						
						
						
						
						

def gen_model(seq_record):
	if(args.multimodel):
		return MultiMarkov.from_string(str(seq_record.seq.upper()), args.order, args.window, maxdiff = 0.001)
	else:
		return GrowingMarkovChain.from_string(str(seq_record.seq.upper()), args.order)	
	
scores = [];	
for seq_record in SeqIO.parse(args.path, "fasta"):
	string = str(seq_record.seq.upper())
	n=1000
	gm = GrowingMarkovChain.from_string(string[:n], args.order)	
	sm = SlidingMarkovChain.from_string(string[n:n*2], args.order)	
	
	for i,(ls, rs) in enumerate(zip(string[n:n*13], string[n*2:n*14])):
		diff = markov_difference(gm, sm)
		scores.append((i, diff));
		gm.grow(ls);
		sm.slide(rs);
		
print max(scores, key=lambda x: x[1])
print scores[8900: 9100]
print string[9500: 10500]
sys.exit();	
	
	
pool = Pool(processes = 2)
res = pool.map(gen_model, SeqIO.parse(args.path, "fasta"))


for mm in res:
	for k, v in mm.emsupport.items():
		print "%s\t%d" % (k, v)
	#for s in mm.generate_string(chunk_size=args.line_length):
		#pass;
	
	#print ">%s" % "_".join(["random", seq_record.id]);
	#for chunk in split2chunks(seq, args.line_length):
		#print chunk;	
	
				
