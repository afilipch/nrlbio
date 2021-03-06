#! /usr/bin/python
'''provides parsing function for stitched maf files and more functions to ananlyse conservation'''
#import miRNA;
import argparse;
import sequencetools;
import parsing;
import os, sys;
from collections import *;

class Maf(object):
	'''Object represents one entry of Maf file
	
	Attributes:
	id str: id of the region(correspondes to id in initial bed file)
	ref_specie str: name of genome assembly(ce6, hg19) for reference specie. Reference specie in alignment does not have gaps
	chr str: chromosome of the region
	strand str: strand of the region
	start int: start of the region
	end int: end of the region
	alignment collections.OrderedDict: represents alignment stack. Key is specie name, Value is sequence from that specie in alignment
	refseq str: sequence of reference specie int alignment block
	'''
	def __init__(self, header, alignment, refseq):
		'''creates a Maf object
		
		header string: main header of an entry in maf file
		alignment collections.OrderedDict: represents alignment stack. Key is specie name, Value is sequence from that specie in alignment
		'''
		self.alignment = alignment;
		self.refseq = refseq;
		
		fl, se, self.id = header.split(":");
		try:
			self.start, self.end = [int(x) for x in se.split(",")[:2]]
		except ValueError:	
			self.start, self.end = [int(x) for x in se.split("-")]##Mackowiak!!!
		sr, chrr = fl.split(".")
		self.ref_specie = sr[1:];
		self.chr, strandr = chrr.split("(")
		self.strand = strandr[:-1];
		
	def __str__(self):
		s = "region ID:%s\nreference_specie:%s\nchromosome:%s\nstrand:%s\nstart:%s\nend:%d\nalignment\n%s%s\n" % (self.id, self.ref_specie, self.chr, self.strand, self.start, self.end, self.ref_specie.ljust(20), self.refseq)
		for k,v in self.alignment.iteritems():
			s+= "%s%s\n" % (k.ljust(20), v)
		return s;

#def read_maf(path, aligned_species=None):
	#'''reads maf stitched file taking into account interactions id. As an output gives dictionary (keys are names of sequence of interest (target seq in interaction human09567). Values are dicts: keys are specie names, values are sequences)
	
	#path string: path to maf stitched file
	#aligned_species list: names of aligned species(genome suystems: ce6, mm9, etc.) to be considered in further analysis. If not given takes all species
	
	#Return dict: Keys are names of sequence of interest (target seq in interaction human09567). Values are dicts: keys are specie names, values are sequences.
	#'''
	#o = defaultdict(dict)
	#f = open(path)
	#for l in f:
		#l = l.strip();  
		#if(l and l[0] == ">"):
			#a = l.split(":");
			#if(len(a) == 3):
				#specie_raw, stub, iid = a
				#specie = specie_raw[1:].split(".")[0];
			#else:
				#specie = l[1:];
		#elif(l and (aligned_species == None or specie in aligned_species)):
			#o[iid][specie] = l.upper();
		#else:
			#pass
	#f.close();
	#sys.stderr.write("%d blocks were read from\t%s\n" % (len(o), path))
	#return o;
	
	
#def read_maf_control(path, aligned_species = None):
	#'''reads maf stitched file. As an output gives dictionary (keys are sequence names (gene names, CLIP cluster names, etc,). Values are dicts: keys are specie names, values are sequences)
	
	#path string: path to maf stitched file
	#aligned_species list: names of aligned species(genome suystems: ce6, mm9, etc.) to be considered in further analysis. If not given takes all species
	
	#Return dict: Keys are sequence names (gene names, CLIP cluster names, etc,). Values are dicts: keys are specie names, values are sequences	
	#'''	
	#maf_control = defaultdict(dict);
	#f = open(path)
	#for l in f:
		#l = l.strip()  
		#if(l and l[0] == ">"):
			#a = l.split(":");
			#if(len(a) > 2):
				#specie_raw, stub = a[:2]
				#gene = ":".join((a[2:]))
				#specie = specie_raw[1:].split(".")[0];
			#else:
				#specie = l[1:];      
		#elif(l and (aligned_species == None or specie in aligned_species)):
			#maf_control[gene][specie] = l.upper();
		#else:
			#pass
	#f.close();
	#return maf_control;	
	

		
def conservation_simple_perfect(match, seqs, pos):
	'''counts how many times given sequence is appeared in alignment stack at the same exact position
	
	match string: sequence to find in each line of alignment (miRNA seed, promoter region, regulatory site and etc.)
	seqs dict: represents alignment stack. Key is specie name, Value is sequence from that specie in alignment
	pos int: position at which match was found in the sequence of reference specie
	
	Return int: how many times given sequence is appeared in alignment stack at the same exact position
	'''
	m = 0;
	end = pos + len(match);
	for k, seq in seqs.iteritems():
		if(match == seq[pos:end]):
			m += 1;
	return	m		
		
	
def conservation_simple_mismatch(match, seqs, pos):
	'''counts how many times given sequence is appeared in alignment stack at the same exact position with one mismatch at any position
	
	match string: sequence to find in each line of alignment (miRNA seed, promoter region, regulatory site and etc.)
	seqs dict: represents alignment stack. Key is specie name, Value is sequence from that specie in alignment
	pos int: position at which match was found in the sequence of reference specie
	
	Return int: how many times given sequence is appeared in alignment stack at the same exact position with one mismatch at any position
	'''	
	m = 0;
	end = pos + len(match);
	variants = sequencetools.diverge_with_1mm(seq, include = True)	
	for k, seq in seqs.iteritems():
		if(seq[pos:end] in variants):
			m += 1;
	return	m;

	
def conservation_perfect(maf, match):
	'''Counts conservation of match in a region represented by Maf object. Requires perfect conservation (same sequence, same position)
	
	maf Maf: Maf object represents region to look for conservation
	match str: sequence to check for conservation (miRNA seed, promoter region, regulatory site and etc.)
	
	Return dict: Keys are positions of match found in reference specie sequence, Values are integers (number of species conserved for a match found in reference specie sequence)	
	'''
	r = {};
	positions = sequencetools.multifind(maf.refseq, match, overlap = False)
	for p in positions:
		r[p] = conservation_simple_perfect(match, maf.alignment, p)
	return r;	
	
	
	
	
	
#def conservation_mismatch(seed, seqs, ref_specie):
	#refseq = seqs[ref_specie];
	#if(seed.match in refseq):
		#return None;
	#for i,level in enumerate(seed.mismatched):
		#for j, mm in enumerate(level):
			#pos = refseq.find(mm.seq)
			#if(pos>-1):
				#r = conservation_simple_mismatch(mm, seed.match, seqs, pos)
				#seed.ass_mismatched[i][j].append(r);	
	#return None;	
	
	
#def conservation_gu_match(seed, seqs, ref_specie):
	#'''function is changing seed object internal property. to be documented'''
	#refseq = seqs[ref_specie];
	#if(seed.match in refseq):
		#pos = refseq.find(seed.match)
		#r = conservation_simple(seed.match, seqs, pos)
		#seed.ass_match.append(r);	
		#return None;	
	#for i,level in enumerate(seed.mismatched):
		#for j, mm in enumerate(level):
			##sys.stderr.write("%s\n" % str(mm.mismatches[0]));
			#if((mm.mismatches[0].fr == "C" and mm.mismatches[0].to == "T") or (mm.mismatches[0].fr == "A" and mm.mismatches[0].to == "G")):
				##sys.stderr.write("bu\n")
				#pos = refseq.find(mm.seq)
				#if(pos>-1):
					#r = conservation_simple_mismatch(mm.seq, seed.match, seqs, pos)
					#seed.ass_match.append(r);	
	#return None;
	
#def conservation_gu_mismatch(seed, seqs, ref_specie):
	#'''function is changing seed object internal property. to be documented'''
	#refseq = seqs[ref_specie];
	#if(seed.match in refseq):
		#return None;
	#for i,level in enumerate(seed.mismatched):
		#for j, mm in enumerate(level):
			#if((mm.mismatches[0].fr == "C" and mm.mismatches[0].to == "T") or (mm.mismatches[0].fr == "A" and mm.mismatches[0].to == "G")):
				#continue;
			#else:		
				#pos = refseq.find(mm.seq)
				#if(pos>-1):
					#r = conservation_simple_mismatch(mm.seq, seed.match, seqs, pos)
					#seed.ass_mismatched[i][j].append(r);
	#return None;	
	
	
  
#def count_conservation_interactions(maf_dict, iid2match, threshold, ref_specie):
	#'''counts conservation along maf_dict aligned species for sequences of interests
	
	#maf_dict dict: Key is id of sequence, Value is dict of aligned sequence (Key is specie name, Value is aligned sequence)
	#iid2match dict: dictionary which connects ids of sequences with actual sequences(miRNA seeds, etc.) to look at conservation. Applicable for analyses of interactions
	#threshold int: min number of species where sequences conserved to consider it as conserved
	#ref_specie str: name of reference specie(for example if one wants to find conservation of human binding site "hg19" should be provided);
	
	#Return tuple: 1st elelment is a dict(Key: sequence of match(miRNA seed match, etc.), Value: number of cases where match appeared to be conserved). 2nd elelment is a dict(Key: sequence of match(miRNA seed match, etc.), Value: number of matches on sequences of reference specie)
	#''' 
	#counts = defaultdict(int);# list of conservation scores
	#total = defaultdict(int);
	#for iid, d in maf_dict.iteritems():
		#match = iid2match[iid]
		#pos = d[ref_specie].find(match)
		#if(pos > -1):
			#nm = conservation_simple(match, d, pos)
			#if(nm >= threshold):
				#counts[match] += 1;
			#total[match] += 1;
	#return counts, total 
  
  
#def count_conservation(maf_dict, matches, threshold):
	#'''counts conservation along maf_dict aligned species for sequences of interests
	
	#maf_dict dict: Key is id of sequence, Value is dict of aligned sequence (Key is specie name, Value is aligned sequence)
	#matches list: list of sequences to check for conservation in each maf_dict block
	#threshold int: min number of species where sequences conserved to consider it as conserved
	#ref_specie str: name of reference specie(for example if one wants to find conservation of human binding site "hg19" should be provided);
	
	#Return tuple: 1st elelment is a dict(Key: sequence of match(miRNA seed match, etc.), Value: number of cases where match appeared to be conserved). 2nd elelment is a dict(Key: sequence of match(miRNA seed match, etc.), Value: number of matches on sequences of reference specie)
	#''' 	
	#counts = defaultdict(int);# list of conservation scores
	#total = defaultdict(int);
	#for g,d in maf_dict.iteritems():
		#for match in matches:      
			#positons = multifind(d[my_specie], match, overlap = False)
			#for pos in positons:
				#nm = conservation_simple(match, d, pos)
				#if(nm >= threshold):
					#counts[match] += 1;
				#total[match] += 1;
	#return counts, total; 	
	
	
	


