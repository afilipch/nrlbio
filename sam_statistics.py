# /usr/bin/python
'''collections of classes and functions to produce various statistics of mapping'''

import os, sys, re;
import Bio;
import pysam;

'''pattern to get matched, mismatched and del/ins stretches to a reference from MD field in sam/bam file file'''
MD_pattern = '([0-9]+)([A-Z]+)*(\^[A-Z]+)*';
MD_compiled_pattern = re.compile(MD_pattern);

def intermediate_alignment(ar):
	'''function to get all mismatches/deletions at defined positions of the aligned read based on the aligned sequence and MD field in sam/bam file. The output is used further to define mismatches/deletions/insertions or/and alignment

	ar pysam.AlignedRead: aligned read to be analyzed

	Return tuple. 1st element is mismatch dictionary (Key: position of mismatch in query, Value: corresponding nucleotide in reference). 2nd element is listiterator (Element: Deleted nucleotides in an order they appear in reference)
	'''	
	mdlist = MD_compiled_pattern.findall(ar.opt('MD'));
	mmdict = {};
	dellist = [];
	p = 0;
	for l in mdlist:
		p += int(l[0]);
		for n in l[1]:
			mmdict[p] = n;
			p+=1;			
		for n in l[2][1:]:
			dellist.append(n);		
		
	deliter = iter(dellist);
	return mmdict, deliter

def get_alignment(ar):	
	'''function to get actual alignment for given pysam.AlignedRead
	
	ar pysam.AlignedRead: aligned read to be analyzed

	Return list of tuples. 1st element in each tuple is nucleotide/gap (string/None) in reference. 2st element in each tuple is nucleotide/gap (string/None) in query.
	'''	
	alignment = [];
	mmdict, deliter = intermediate_alignment(ar);
	ins_adjust = 0;	
	for i, j in ar.aligned_pairs:
		if(i != None):
			if(j != None):
				if(i - ins_adjust in mmdict):
					apair = (mmdict[i-ins_adjust], ar.query[i]); 
				else:
					apair = (ar.query[i], ar.query[i]);
			else:
				ins_adjust += 1;
				apair = (None, ar.query[i])
		else:
			apair = (next(deliter), None)
		alignment.append(apair);	
			
	return alignment;
	
def get_conversions(ar):
	'''function to get all mismatches/deletions/insertions for given pysam.AlignedRead
	
	ar pysam.AlignedRead: aligned read to be analyzed

	Return list of tuples. 1st element in each tuple is nucleotide/gap (string/None) in reference. 2st element in each tuple is nucleotide/gap (string/None) in query.
	'''	
	conversions = [];
	mmdict, deliter = intermediate_alignment(ar);
	ins_adjust = 0;	
	for i, j in ar.aligned_pairs:
		if(i != None):
			if(j != None):
				if(i - ins_adjust in mmdict):
					apair = (mmdict[i-ins_adjust], ar.query[i]); 
				else:
					continue;
			else:
				ins_adjust += 1;
				apair = (None, ar.query[i])
		else:
			apair = (next(deliter), None)
		conversions.append(apair);	
			
	return conversions;	
	
	

