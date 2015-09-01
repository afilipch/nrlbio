#! /usr/bin/python
import sys;
import parsing;
import re;
import os;
import time;
import csv;
import copy;
import argparse;
import random
import multiprocessing;
from collections import *;
from Bio import SeqIO;
from generators import *;
from Bio import pairwise2;

parser = argparse.ArgumentParser(description='search for highly enriched mers in the reads, outputs results to  STDOUT');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to fastq file of reads");
parser.add_argument('-b', '--backward', nargs = '?', default = False, const = True, type = str, help = "search for 5'adapter");
parser.add_argument('-l', '--length', nargs = '?', default = 10, type = int, help = "length of mer");
parser.add_argument('-m', '--mode', nargs = '?', default = "all", type = str, help = "search for all kmers, from start, from end; can be all|start|end");
parser.add_argument('-n', '--steps', nargs = '?', default = 80, const = 1, type = int, help = "number of steps from start or end position");
parser.add_argument('-a', '--adapters', nargs = '?', required = True, help = "output alignments of the most enriched mers to the adapters provided in fasta format");
parser.add_argument('--most_common', nargs = '?', default = 20, type = int, help = "number of top found subsequences to be shown");
args = parser.parse_args();

class Anchor(object):
	"""represent anchor 5,6,7 etc mer that further can be extended"""
	def __init__(self, anchor, tail, position):
		"""anchor: sequence of nmer to be anchor extracted form ref.seq. in the read sequence, tail: lef or right tail of ref.seq. adjacent to the anchor, postition: position of nmer inside ref.seq. start of nmer for forward end of nmer for backward, self.sequences: defaultdict(int) counting presence of individual sequences (anchor + part of tail) in reads"""
		self.anchor = anchor;
		self.tail = tail;
		self.position = position;
		#self.sequences = defaultdict(int);
		
	def extend(self, read_sequence):
		"""try to extend anchor in forward direction. as input takes right part of read sequence adjacent to anchor. as ouput increment the value in self.sequences corresponding to longest extension found""" 
		a = "";
		for i in range(min(len(read_sequence), len(self.tail))):
			if(read_sequence[i] == self.tail[i]):
				a += read_sequence[i];
			else:
				#self.sequences[self.anchor + a] += 1;
				return self.anchor + a, i;
		#self.sequences[self.anchor + a] += 1;
		return	self.anchor + a, min(len(read_sequence), len(self.tail));

class RefSeq(object):
	"""represent reference sequence, allows to generate anchors of given sequence"""
	def __init__(self, sequence, name):
		self.sequence = sequence;
		self.name = name;
		self.anchors = {};
		self.sequences = defaultdict(list);
		
	def generate_anchors(self, length):
		for i in range(len(self.sequence) - length + 1):
			if(self.sequence[i:i+length] in self.anchors):
				continue;
			else:	
				self.anchors[self.sequence[i:i+length]] = Anchor(self.sequence[i:i+length], self.sequence[i+length:], i)	
				
	def get_sequences(self):
		"""produces all possible subsequences of reference sequence"""
		for anchor in self.anchors.values():
			a = anchor.anchor;
			self.sequences[anchor.position].append(a)
			for char in anchor.tail:
				a += char;
				self.sequences[anchor.position].append(a);
				


def search(refseq, sequence, length, backward):
	"""search for refseq nonitersected subsequences inside the sequence(read)"""
	if(backward):
		sequence = sequence[::-1]
	pos = 0;
	pieces = [];
	while (pos < len(sequence) - length + 1):
		mer = sequence[pos:pos+length]
		anchor = refseq.anchors.get(mer, None);
		if(anchor):
			piece, tpos = anchor.extend(sequence[pos+length:]);
			pieces.append(piece);
			pos = pos + tpos + length;
		else:
			pos += 1;
	return pieces;		
	
	
def search_inlist(arr):
	"""main function outputs Counter objects with adapter subsequences as keys, and number of reads caontaining them, nonintersected pieces in adapter sequence are counted separately"""
	refseq_list, seq_list, length, backward = arr
	ans = {};
	#print len(seq_list)
	
	for refseq in refseq_list:
		ans[refseq.name] = Counter();
		for seq in seq_list:
			sequence = seq[0];
			
			pos = 0;
			while (pos < len(sequence) - length + 1):
				mer = sequence[pos:pos+length]
				anchor = refseq.anchors.get(mer, None);
				if(anchor):
					piece, tpos = anchor.extend(sequence[pos+length:]);
					ans[refseq.name][piece] += 1;
					pos = pos + tpos + length;
				else:
					pos += 1;			
	return ans;	
		
def represent(refseq, common, backward):
	print "_"*100
	my_common = [];
	if(backward):
		seq = refseq.sequence[::-1];
		for k,v in common:
			my_common.append((k[::-1], v))
	else:
		seq = refseq.sequence;
		my_common = common;
		
	for k,v in my_common:
		if(backward):
			pos = seq.rfind(k)
		else:
			pos = seq.find(k)
			
		print pos*" " + k + "\t" + str(v)		
	print
	print seq + "\t" + refseq.name
	print 
	
	
## Workflow	
	
#adapters are read as dictionary from file provided to "--adapters" option  
adapters = parsing.fasta2dict(args.adapters, reverse = args.backward);

#create list of refseq objects on basis of adapters fasta dictionary
refseq_list = []
for k, v in adapters.iteritems():
	refseq = RefSeq(v, k);
	refseq.generate_anchors(args.length)
	refseq.get_sequences()
	refseq_list.append(refseq);


# keys: names of adapters, values:  Counter objects for each adapters: adapter subsequences as keys, and number of reads caontaining them, nonintersected pieces in adapter sequence are counted separately
pieces	= {}
for refseq in refseq_list:
	pieces[refseq.name] = Counter();

	
pool = multiprocessing.Pool(processes = 5)
res = pool.imap(search_inlist, ((refseq_list, seq_list, args.length, args.backward) for seq_list in grouper(generator_fastq(args.path, ["seq"], reverse = args.backward), 10000)) )      
for temp in res:
	for k, v in temp.iteritems():
		pieces[k] += v;

		
for refseq in refseq_list:  
	represent(refseq, pieces[refseq.name].most_common(args.most_common), args.backward);	






















