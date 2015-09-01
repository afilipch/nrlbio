import parsing;
import sys;
import os;
import time;
from itertools import *;
from collections import *;

def get_mode(s):
	if(s == "hc"):
		return {};
	elif(s == "pc"):	
		return {"T": ["C"]}
	elif(s == "exp"):			
		exp_mode = defaultdict(list);
		for n1 in "ACTG":
			for n2 in ["A", "T", "C", "G", ""]:
				if(n1 != n2):
					exp_mode[n1].append(n2);
		return exp_mode
	else:
		raise Exception("invalid input value")
			

def mutate_piece(piece, mode):	
	'''from initial sequence creates many edited with mutations listed in mode argument'''
	my_pieces = set();
	my_pieces.add(piece);
	for i in range(len(piece)):
		if(piece[i] in mode):
			for n in mode[piece[i]]:
				my_pieces.add(piece[:i] + n + piece[i+1:]);
	return my_pieces;  

 
def get_pieces(ref_dict, first, length, mode):
	'''creates edited pieces form given sequences'''
	pieces = defaultdict(set);
	for rid, seq in ref_dict.iteritems():
		my_pieces = set();
		for i in range(min(first, len(seq) - length + 1)):
			piece = seq[i: i + length];
			my_pieces.update(mutate_piece(piece, mode));
			
		for pie in my_pieces:
			pieces[pie].add(rid); 
	sys.stderr.write("total pieces generated %d\n" % len(pieces))		
	return pieces;
	

    
    
    
    
