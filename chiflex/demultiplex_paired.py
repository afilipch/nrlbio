#! /usr/lib/python
'''Demultiplexing of raw mapping reasults for paired-end sequncing'''
import argparse
import os;
import sys
from collections import defaultdict
from itertools import product, combinations

import pysam;

from nrlbio.statistics import sam as sam_statistics
from nrlbio.samlib import ArWrapper, demultiplex_read_hits




parser = argparse.ArgumentParser(description='Assignes to each read if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', default = "sam", type = str, help = "path to the output folder");
parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
parser.add_argument('-s', '--score', nargs = '?', default = 'as', choices = ['as', 'as_qstart', 'as_qstart_entropy', 'as_qstart_rstart', 'as_qstart_rstart_entropy'], type = str, help = "score function for hits");
parser.add_argument('--maxdistance', nargs = '?', default = 1000, type = int, help = "Max distance allowed for paired-end hits to be considered as concordant alignment");
parser.add_argument('--maxgap', nargs = '?', default = 0, type = int, help = "Max gap between hits on the same read allowed to consider this read(one mate) to be chimeric");
parser.add_argument('--maxoverlap', nargs = '?', default = 6, type = int, help = "Max overlap between hits on the same read allowed to consider this read(one mate) to be chimeric");
parser.add_argument('--concreward', nargs = '?', default = 3, type = float, help = "Multiplier for the score for concordant mappings");
args = parser.parse_args();

# parse [score] argument
score2function = {'as': 'as_score', 'as_qstart': 'as_qstart_score', 'as_qstart_entropy': 'as_qstart_entropy_score', 'as_qstart_rstart': 'as_qstart_rstart_score', 'as_qstart_rstart_entropy': 'as_qstart_rstart_entropy_score'}
exec("from nrlbio.samlib import %s as key_score" % score2function[args.score]);


def doublet_score(arws):
	return sum([x.score for x in arws])

def triplet_score(ch1, ch2, arw):
	return ch1.score + ch2.score + arw.score;


def clear_dict(d):
	for v in d.values():
		v[:]=[];
		
def mate2chimera(arws, maxoverlap, maxgap):
	distance = max([x.qstart for x in arws]) - min([x.qend for x in arws])
	if(-1*maxoverlap<=distance<=maxgap):
		return tuple(sorted(arws, key = lambda x: x.qstart))
	else:
		return None
	
def triplet2chimera(matechimeras, arws, maxdistance, concreward):
	triplets = [];
	for (ch1,ch2), arw in product(matechimeras, arws):
		if(ch2.rname == arw.rname):
			distance = max(ch2.aligned_read.reference_start, arw.aligned_read.reference_start) - min(ch2.aligned_read.reference_end, arw.aligned_read.reference_end)
			if(distance<=maxdistance):
				triplets.append(((ch1, ch2, arw), triplet_score(ch1, ch2, arw)*concreward));
	return triplets;
	
	
def doublet2chimera(arws, maxdistance, concreward):
	if(arws[0].rname == arws[1].rname):
		distance = max([x.aligned_read.reference_start for x in arws]) - min([x.aligned_read.reference_end for x in arws])
		if(distance<=maxdistance):
			return (arws, doublet_score(arws)*concreward);
		else:
			return (arws, doublet_score(arws));
	else:
		return (arws, doublet_score(arws));
	

def demultiplex(paired_segments, maxdistance, maxoverlap, maxgap, concreward):
	doublets = [];
	for arws in product(paired_segments['1'], paired_segments['2']):
		doublets.append(doublet2chimera(arws, maxdistance, concreward))
	if(doublets):
		best_doublet, doublet_score = max(doublets, key = lambda x: x[1])
	else:
		best_doublet, doublet_score = None, -1;
		
	#if(best_doublet and not any([x.control for x in best_doublet])):
		#l = [x.rname.split("|")[-1] for  x in best_doublet]
		#if(l and l[0]!=l[1]):
			#l.sort();
			#print "\t".join(l)
			
	matechimeras = defaultdict(list)
	for k, v in paired_segments.items():
		for arws in combinations(v, 2):
			mc = mate2chimera(arws, maxoverlap, maxgap)
			if(mc):
				matechimeras[k].append(mc);
				l = sorted([x.rname.split("|")[-1] for  x in mc])
				if(l[0]!=l[1]):
					print "\t".join(sorted([x.rname.split("|")[-1] for  x in mc]))
				
	#112 triplet case
	triplets = triplet2chimera(matechimeras['1'], paired_segments['2'], maxdistance, concreward)
	triplets += triplet2chimera(matechimeras['2'], paired_segments['1'], maxdistance, concreward)
	if(triplets):
		best_triplet, triplet_score = max(triplets, key = lambda x: x[1])
	else:
		best_triplet, triplet_score = None, -1;
			
				
	#if(best_triplet and not any([x.control for x in best_triplet])):
		#l = list(set([x.rname.split("|")[-1] for  x in best_triplet]))
		#if(l and len(l)==2):
			#l.sort();
			#print "\t".join(l)		
	
#######################################test functions	
	
def just_pairs(paired_segments):
	'''function for first step comparison'''
	#print '*'*140
	for p in product(*paired_segments.values()):
		l = [x.rname.split("|")[-1] for  x in p]
		if(l and l[0]!=l[1]):
			l.sort();
			print "\t".join(l)	
	
	
def check_reverse(paired_segments, maxdistance, maxoverlap, maxgap, concreward):
	for arw in paired_segments['1']:
		print arw.aligned_read.is_reverse
	
demultiplex = check_reverse	
	


samfile = pysam.Samfile(args.path);
paired_segments = defaultdict(list);
current_name = '';

for segment in samfile.fetch(until_eof=True):
	if(not segment.is_unmapped):
		rname = samfile.getrname(segment.tid)
		arw = ArWrapper(segment, rname, score_function=key_score, add_nr_tag=False)
		basename, mate = segment.query_name.split(':');
		if(current_name != basename):
			demultiplex(paired_segments, args.maxdistance, args.maxoverlap, args.maxgap, args.concreward);
			clear_dict(paired_segments);
			current_name = basename;
		paired_segments[mate].append(arw);
#else:
	#demultiplex(paired_segments, maxdistance, maxoverlap, maxdgap, concreward);
	
		
		
		
		
		