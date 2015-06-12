#! /usr/lib/python
'''Generates decoy fasta file on basis of given one'''

import argparse
import os
import sys
import collections
from multiprocessing import Pool

from Bio import SeqIO

from nrlbio.sequencetools import split2chunks;
from nrlbio.HMM import MarkovChain, MultiMarkov, SlidingMarkovChain, GrowingMarkovChain, markov_difference

parser = argparse.ArgumentParser(description='Generates decoy fasta file on basis of given one');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file");
parser.add_argument('-m', '--multimodel', nargs = '?', default = False, const = True, type = bool, help = "Use set of Markov models to generate decoy sequence");
parser.add_argument('-o', '--order', nargs = '?', default = 1, type = int, help = "Markov model order. Correspondes to shuffling of \'order\'-nucleotides");
parser.add_argument('-w', '--window', nargs = '?', default = 400, type = int, help = "Can be used in multimodel regime. Each fasta sequence will be split into sequences of length \'window\', for each of them separate Markov Model will be generated");
parser.add_argument('-ll', '--line_length', nargs = '?', default = 100, type = int, help = "Max length of line in output fasta file. Makes it more readable");
args = parser.parse_args();


#def find_edge(string, growing, sliding, slide_size, jump_diff, max_diff, curpos):
	#jump = slide_size/4
	#max_look = 100;
	#maxscore = 0;
	#lookforward = 0;
	
	
	
	
	#for i, (ls, rs) in enumerate(zip(string, string[slide_size:])):
		#diff = markov_difference(growing, sliding);
		
		#if(lookforward==max_look and maxscore>max_diff):
			#growing = growing.shrink(string[i-lookforward-growing.order*2:i])
			#ngrowing = GrowingMarkovChain.from_string(string[i-lookforward: i+slide_size-lookforward], args.order)
			#nsliding = SlidingMarkovChain.from_string(string[i+slide_size-lookforward: i+slide_size*2-lookforward], args.order)
			
			#sys.stderr.write("switch to a new model happens at position (%d), cause the difference (%1.5f) is more than max difference %1.5f and lookforward is more or equal than %d\n" % 
			#(i+curpos-lookforward, maxscore, max_diff, lookforward))
			#return growing, ngrowing, nsliding, i-lookforward+slide_size;
			
		#elif(diff<jump_diff or lookforward>max_look*5):
			#growing.grow_long(string[i: i+jump])
			#sliding = SlidingMarkovChain.from_string(string[i+jump: i+slide_size+jump], sliding.order)
			#sys.stderr.write("jump of length (%d) happens at position (%d), cause the difference (%1.5f) is less than jump difference %1.5f\n" % (jump, i+curpos, diff, jump_diff))
			#return None, growing, sliding, i+jump;
			
		#else:
			#growing.grow(ls);
			#sliding.slide(rs);	
			#if(diff>maxscore):
				#lookforward = 0;
				#maxscore = diff;
			#else:
				#lookforward +=1;
	#else:
		#growing.add(sliding)
		#return growing, None, None, 0
		
		
def seqrecord2models(seqrecord, step_size=50000, slide_size=1200):
	models = []
	steps = len(seqrecord)/step_size;
	
	sz=1200
	curpos=sz
	ngrowing = GrowingMarkovChain.from_string(string[:sz], args.order);
	nsliding = SlidingMarkovChain.from_string(string[sz:sz*2], args.order);

	while(ngrowing):
		growing, ngrowing, nsliding, pos = find_edge(string[curpos:], ngrowing, nsliding, slide_size=sz, jump_diff= 0.009, max_diff=0.01, curpos=curpos);
		curpos += pos;
		if(growing):
			models.append(growing);		
		

models = [];
step = 50000
for seq_record in SeqIO.parse(args.path, "fasta"):
	k=0;
	lseq = len(seq_record);
	print lseq;
	while(k+step<lseq):
		print k
		string2models(str(seq_record[k:k+step].seq.upper()));
		k+=step
			
print len(models)		
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
	
				
