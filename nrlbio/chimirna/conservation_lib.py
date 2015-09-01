#! /usr/bin/python
'''provides parsing function for stitched maf files and more functions to ananlyse conservation'''
import mir_lib;
import argparse;
import parsing;
import os, sys;
from collections import *;

def read_maf(path, aligned_species=None):
	'''reads maf stitched file taking into account interactions id. As an output gives dictionary (keys are names of sequence of interest (target seq in interaction human09567). Values are dicts: keys are specie names, values are sequences)
	
	path string: path to maf stitched file
	aligned_species list: names of aligned species(genome suystems: ce6, mm9, etc.) to be considered in further analysis. If not given takes all species
	
	Return dict: Keys are names of sequence of interest (target seq in interaction human09567). Values are dicts: keys are specie names, values are sequences.
	'''
	o = defaultdict(dict)
	f = open(path)
	for l in f:
		l = l.strip();  
		if(l and l[0] == ">"):
			a = l.split(":");
			if(len(a) == 3):
				specie_raw, stub, iid = a
				specie = specie_raw[1:].split(".")[0];
			else:
				specie = l[1:];
		elif(l and (aligned_species == None or specie in aligned_species)):
			o[iid][specie] = l.upper();
		else:
			pass
	f.close();
	sys.stderr.write("%d blocks were read from\t%s\n" % (len(o), path))
	return o;
	
	
def read_maf_control(path, aligned_species = None):
	'''reads maf stitched file. As an output gives dictionary (keys are sequence names (gene names, CLIP cluster names, etc,). Values are dicts: keys are specie names, values are sequences)
	
	path string: path to maf stitched file
	aligned_species list: names of aligned species(genome suystems: ce6, mm9, etc.) to be considered in further analysis. If not given takes all species
	
	Return dict: Keys are sequence names (gene names, CLIP cluster names, etc,). Values are dicts: keys are specie names, values are sequences	
	'''	
	maf_control = defaultdict(dict);
	f = open(path)
	for l in f:
		l = l.strip()  
		if(l and l[0] == ">"):
			a = l.split(":");
			if(len(a) > 2):
				specie_raw, stub = a[:2]
				gene = ":".join((a[2:]))
				specie = specie_raw[1:].split(".")[0];
			else:
				specie = l[1:];      
		elif(l and (aligned_species == None or specie in aligned_species)):
			maf_control[gene][specie] = l.upper();
		else:
			pass
	f.close();
	return maf_control;	
	
def conservation(seed, dictionary, pos, matches):
	m = 0;
	for key, seq in dictionary.iteritems():
		if(seed == seq[pos:pos+6]):
			m += 1;
	if(m <= matches):
		return 0;
	else:
		return 1;
		
def conservation_simple(match, seqs, pos):
	m = 0;
	end = pos + len(match);
	for k, seq in seqs.iteritems():
		if(match == seq[pos:end]):
			m += 1;
	return	m	
	
def conservation_perfect(seed, seqs, ref_specie):
	refseq = seqs[ref_specie];
	pos = refseq.find(seed.match)
	if(pos>-1):
		r = conservation_simple(seed.match, seqs, pos)
		seed.ass_match.append(r);	
	return None;	
		
	
def conservation_simple_mismatch(mm, match, seqs, pos):
	m = 0;
	end = pos + len(match);
	variants = [];
	for mposition in range(len(match)):
		for n in "ACTG":
			variants.append(match[:mposition] + n + match[mposition+1:])	
	#for n in "ACTG":
		#variants.append(match[:mm.mismatches[0].position] + n + match[mm.mismatches[0].position+1:])
	for k, seq in seqs.iteritems():
		if(seq[pos:end] in variants):
			m += 1;
	return	m	
	
def conservation_mismatch(seed, seqs, ref_specie):
	refseq = seqs[ref_specie];
	if(seed.match in refseq):
		return None;
	for i,level in enumerate(seed.mismatched):
		for j, mm in enumerate(level):
			pos = refseq.find(mm.seq)
			if(pos>-1):
				r = conservation_simple_mismatch(mm, seed.match, seqs, pos)
				seed.ass_mismatched[i][j].append(r);	
	return None;	
	
	
def conservation_gu_match(seed, seqs, ref_specie):
	refseq = seqs[ref_specie];
	if(seed.match in refseq):
		pos = refseq.find(seed.match)
		r = conservation_simple(seed.match, seqs, pos)
		seed.ass_match.append(r);	
		return None;	
	for i,level in enumerate(seed.mismatched):
		for j, mm in enumerate(level):
			#sys.stderr.write("%s\n" % str(mm.mismatches[0]));
			if((mm.mismatches[0].fr == "C" and mm.mismatches[0].to == "T") or (mm.mismatches[0].fr == "A" and mm.mismatches[0].to == "G")):
				#sys.stderr.write("bu\n")
				pos = refseq.find(mm.seq)
				if(pos>-1):
					r = conservation_simple_mismatch(mm.seq, seed.match, seqs, pos)
					seed.ass_match.append(r);	
	return None;
	
def conservation_gu_mismatch(seed, seqs, ref_specie):
	refseq = seqs[ref_specie];
	if(seed.match in refseq):
		return None;
	for i,level in enumerate(seed.mismatched):
		for j, mm in enumerate(level):
			if((mm.mismatches[0].fr == "C" and mm.mismatches[0].to == "T") or (mm.mismatches[0].fr == "A" and mm.mismatches[0].to == "G")):
				continue;
			else:		
				pos = refseq.find(mm.seq)
				if(pos>-1):
					r = conservation_simple_mismatch(mm.seq, seed.match, seqs, pos)
					seed.ass_mismatched[i][j].append(r);
	return None;	
	
	
  
def count_conservation(maf_dict, iid2match, threshold, ref_specie):
	'''counts conservation along maf_dict aligned species for sequences of interests
	maf_dict Dict: dictionary (keys: ids of sequences, values: dicts of aligned sequences)
	iid2match Dict: dictionary which connects ids of sequences with actual sequences to look at conservation
	threshold Int: min number of species where sequences conserved to consider it as conserved
	ref_specie String: name of reference specie(for example if one wants to find conservation of human binding site "hg19" should be provided);
	''' 
	counts = defaultdict(int);# list of conservation scores
	total = defaultdict(int);
	for iid, d in maf_dict.iteritems():
		match = iid2match[iid]
		pos = d[ref_specie].find(match)
		if(pos > -1):
			counts[match] += conservation(match, d, pos, threshold)
			total[match] += 1;
	return counts, total 
  
  
def count_control(maf_control, matchs, threshold):
	'''counts conservation along maf_dict aligned species for sequences of interests
	maf_dict Dict: dictionary (keys: ids of sequences, values: dicts of aligned sequences)
	matches Lict: list of sequences to check for conservation in each maf_dict block
	threshold Int: min number of species where sequences conserved to consider it as conserved
	ref_specie String: name of reference specie(for example if one wants to find conservation of human binding site "hg19" should be provided);
	''' 	
	counts = defaultdict(int);# list of conservation scores
	total = defaultdict(int);
	for g,d in maf_control.iteritems():
		for match in matchs:      
			positons = parsing.findsubstring(d[my_specie], match)
			for pos in positons:
				counts[match] += conservation(match, d, pos, threshold)
				total[match] += 1;
	return counts, total; 	
	
	
	


