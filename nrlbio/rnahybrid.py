'''Library contains some wrappers/parsers for RNA hybridization tools'''
import commands;
import random;
import sys;
from collections import defaultdict

from nrlbio.sequencetools import shuffle_string


### RNAhybrid wrappers

def _paired_rnahybrid(basepairing):
	'''Extracts pairing information from raw RNAhybrid text output 
		basepairing list: strings representing raw RNAhybrid text output 
	Returns list: list of values 0-> unbound nucleotide,1-> bound nucleotide along miRNA
	'''
	bound = [];
	micro_match = basepairing[2][::-1]; # 5->3 direction
	micro_unmatch = basepairing[3][::-1]; # 5->3 direction
	for i in range(len(micro_match)):
		if(micro_match[i] != " "):
			bound.append(1);
		elif(micro_unmatch[i] != " "):
			bound.append(0);
	return bound;  


def get_rnahybrid(target, mirna, system='3utr_human', arguments = ""):
	"""Calculates hybridization energy of intermolecular interaction between given sequences via RNAhybrid tool. NOTE, RNAhybridis is designed to hybridize miRNA and respective target (not folding themselves). Therefore input sequences are called 'target' and 'mirna' respectively, and basepairing of mirna is also reported.
	
		target str: miRNAs target sequence
		mirna str: miRNA sequence
		system str: background system for RNAhybrid, does not affect current workflow of the function
		arguments str: optional arguments for RNAhybid. For detailed information see RNAhybid documentation
		
	Returns: 
		float: hybridization energy
		list: basepairing, 0-> unbound nucleotide,1-> bound nucleotide along miRNA
	"""
	call =  "RNAhybrid -s %s  \"%s\" \"%s\" -b 1 -c %s " % (system, target, mirna, arguments)
	status, output = commands.getstatusoutput(call)
	e = float(output.split(":")[4]);  
	basepairing = output.split(":")[7:11];  
	p = _paired_rnahybrid(basepairing);
	return e, p;


def get_shuffled_rnahybrid(target, mirna, trials=20, system='3utr_human', arguments = ""):
	"""Calculates background(sequence shuffling) hybridization energy of intermolecular interaction between given sequences via RNAhybrid tool. NOTE, RNAhybridis is designed to hybridize miRNA and respective target (not folding themselves). Therefore input sequences are called 'target' and 'mirna' respectively, and basepairing of mirna is also reported. NOTE: mean energy and hybridiztion pattern of different shuffling trials is returned. 
	
		target str: miRNAs target sequence
		mirna str: miRNA sequence
		trials int: controls how many times sequences will be shuffled. NOTE: mean energy and hybridiztion pattern of different shuffling trials is returned.
		system str: background system for RNAhybrid, does not affect current workflow of the function
		arguments str: optional arguments for RNAhybid. For detailed information see RNAhybid documentation
		
	Returns: 
		float: hybridization energy
		list: basepairing, 0-> unbound nucleotide,1-> bound nucleotide along miRNA
	"""
	
	energies = [];
	patterns = defaultdict(list);
	
	for c in range(trials):
		e, pattern = get_rnahybrid(shuffle_string(target), mirna, system=system, arguments = arguments);
		energies.append(e);
		for c, p in enumerate(pattern):
			patterns[c].append(p);
			
	return sum(energies)/float(trials), [sum(patterns[x])/float(trials) for x in range(len(patterns))]
			


