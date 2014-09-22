'''basic function to deal with nucleotide sequence'''
import random;
import re;
import sys

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
	for p in range(len(match)):
		for n in nset - set(match[p]):
			variants.append(match[:mposition] + n + match[mposition+1:]);
	return variants;	
	
	
def shuffle_string(s):
	'''shuffles given string'''
	l = list(s)
	random.shuffle(l)
	return ''.join(l)
	
	
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
		
		
def entropy(seq):
	sys.stderr.write('Is not implemented yet!\n')
		
		
		
if(__name__ == "__main__"):
	seq = "ANDREI"
	length = 5;
	for chunk in split2chunks(seq, length):
		print chunk;
			