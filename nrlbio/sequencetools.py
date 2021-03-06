'''basic function to deal with nucleotide sequence'''
import random;
import re;
import sys;
from collections import defaultdict, Counter;
from numerictools import dict2entropy
from itertools import combinations, product

RevComplDict = {'G': 'C', 'C': 'G', 'A': 'T', 'T': 'A'}
NUCLEOTIDES = set("ACTG")

def diverge_with_1mm(seq, include = False):
	'''produces all sequences which are 1nt different from the initial one
	
	seq string: nucleotide sequence to be mutated
	include bool: if True initial sequence is included in output list
	
	Return list: all sequences which are 1nt different from the initial one
	'''
	if(include):
		variants = [seq];
	else:	
		variants = [];
	nset = set("ACTG");
	for p in range(len(seq)):
		for n in nset - set(seq[p]):
			variants.append(seq[:p] + n + seq[p+1:]);
	return variants;
	
	
def generate_mismatched_sequence(seq, number_of_mismatches = 1):
	'''generates sequences with the exact number of mismatches to the initial one
	
		seq str: initial sequence to generate mismatched from
		number_of_mismatches int: number of mismatches in generated sequence compare to the initial seq

		Yields str: mismatched sequence
	'''	
	
	for positions in combinations(range(len(seq)), number_of_mismatches):
		
		variants = [];
		for p in positions:
			locker = set(seq[p]);
			variants.append(NUCLEOTIDES-locker)
			
		mm_at_positions = product(*variants) 
		
		for mm in mm_at_positions:
			l = list(seq)
			
			for i, m in enumerate(mm):
				l[positions[i]] = m;
				
			yield "".join(l)
			
			
def introduce_conversions(seq, probabilty):
	'''intoduces random conversions into a given string with defined probability
	
	Returns str, int: mutated sequence, number of conversions
	'''
	nlist = [];
	num = 0;
	for s in seq:
		if(random.random()<probabilty):
			nlist.append(random.choice(list(NUCLEOTIDES-set(s))))
			num +=1;
		else:
			nlist.append(s);
	return "".join(nlist), num


	
	
def shuffle_string(s):
	'''shuffles given string'''
	l = list(s)
	random.shuffle(l)
	return ''.join(l)
	
	
def random_string(length):
	'''Returns randomly generated nucleotide sequences'''
	l = []
	for _ in range(length):
		l.append(random.choice('ACTG'));
	return ''.join(l);	
		
		
def generate_random_strings(length, number):	
	'''Yields randomly generated nucleotide sequences'''
	for _ in range(number):
		yield random_string(length);	
		
		
	
def multifind(string, substring, overlap = False):
	'''finds all occurences of substring in the string
	
	string str: string to find substring in
	substring str: substring to be found inside the string
	overlap bool: if True, function looks for overlapping substring
	
	Return list: list of start positions of substring inside the string. If no matches found, empty list is returned
	'''	
	if(overlap):
		p = '(?=' + substring + ')'
	else:
		p = substring
	return 	[m.start() for m in re.finditer(p, string)]
	
	
def split2chunks(seq, length):
	""" Yield successive length-sized chunks from seq"""
	for i in xrange(0, len(seq), length):
		yield seq[i:i+length]
		

		
def _get_transitions(seq, order):
	if(order):
		transitions = defaultdict(int)
		for i in range(len(seq)-2*order+1):
			transitions[ seq[i:i+order], seq[i+order:i+2*order] ] += 1;
	else: 
		transitions = Counter(seq);
	return transitions
		
def entropy(seq, order=1):
	'''Calculates Shannon entropy of provided sequence on basis of Markov Model transition probabilities
	
		seq str: sequence to calculate entropy of
		order int: order of underlying Markov Model
		
	Returns float: entropy 	of provided sequence
	'''
	return dict2entropy(_get_transitions(seq, order))
	
	
def chunk_entropy(seq, length, step = 1, order = 1):	
	'''Returns minimum Shannon entropy of chunks from provided sequence on basis of Markov Model transition probabilities
	
		seq str: sequence to calculate entropy of
		length int: length of chunk to calculate entropy
		step int: step for the sliding window to generate pieces
		order int: order of underlying Markov Model
		
	Returns float: entropy 	of provided sequence
	'''
	if(len(seq) < length):
		return 0;
		
	min_entropy = dict2entropy(_get_transitions(seq[:length], order=order));
	for i in range(step, len(seq)-length+1, step):
		min_entropy = min(min_entropy, dict2entropy(_get_transitions(seq[i:length+i], order=order)));
	return	min_entropy


def splitstring(s, size):
	'''splits string into chunks of particular length'''
	return [s[x:x+size] for x in range(0, len(s), size)]


	
		
		
if(__name__ == "__main__"):
	#rep1 = 'AAAAAAAATAAA'
	#rep2 = 'ATATATATATAT'
	#seq1 = 'CGTCATCAAGCA'
	#sel1 = 'TCACATGACTAGCGTCATCAAGCAGCATGCGTACACAGTCAGTCAACAGAGCAGATATTCAAATCAGCTAGGCACCATGACGCTATATATGGGGGG'
	#print entropy(rep1)
	#print entropy(rep2)
	#print entropy(rep1+rep2)
	#print entropy(seq1)
	#print entropy(sel1)	
	#print chunk_entropy(sel1, 12, step=1)
	
	#n = 0;
	#for seq in random_string(50, 4000):
		##print seq
		#me = chunk_entropy(seq, 12, step=1);
		#if(me>1.25):
			#n+=1;
	#print n
	
	#for s in generate_mismatched_sequence("AGT", number_of_mismatches = 2):
		#print s;
	
	
	print introduce_conversions('AAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 0.1)
			
			
			
			
			