#! /usr/lib/python
'''Demultiplexing of raw mapping reasults for paired-end sequncing'''
import argparse
import os;
import sys
from collections import defaultdict, namedtuple
from itertools import product, combinations, izip

import pysam;

from nrlbio.statistics import sam as sam_statistics
from nrlbio.samlib import ArWrapper, demultiplex_read_hits
from nrlbio.numerictools import maxes
from nrlbio.pybedtools_extension import construct_gff_interval
from nrlbio.generators import generator_segments




parser = argparse.ArgumentParser(description='Assignes to each read if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "path to sam/bam files. The first file is a mapping of the first mate, the second is for the second mate");
parser.add_argument('-o', '--output', nargs = '?', default = "chimeras", type = str, help = "path to the output folder");
parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for output files, should reflect nature of mapping reference");
parser.add_argument('-s', '--score', nargs = '?', default = 'as', choices = ['as', 'as_qstart', 'as_qstart_entropy', 'as_qstart_rstart', 'as_qstart_rstart_entropy'], type = str, help = "score function for hits");
parser.add_argument('--bestdistance', nargs = '?', default = 8, type = float, help = "minimal distance allowed between the best and the second best hit. If the actual distance is less, than hit will be assigned as nonunique");
parser.add_argument('--maxdistance', nargs = '?', default = 1000, type = int, help = "Max distance allowed for paired-end hits to be considered as concordant alignment");
parser.add_argument('--maxgap', nargs = '?', default = 0, type = int, help = "Max gap between hits on the same read allowed to consider this read(one mate) to be chimeric");
parser.add_argument('--maxoverlap', nargs = '?', default = 6, type = int, help = "Max overlap between hits on the same read allowed to consider this read(one mate) to be chimeric");
parser.add_argument('--concreward', nargs = '?', default = 2, type = float, help = "Multiplier for the score for concordant mappings");
parser.add_argument('--verbose', nargs = '?', default = 0, const = 1000000, type = int, help = "If set, progress of the demultiplexing will be reported");
args = parser.parse_args();

##Set constants
Chimera = namedtuple('Chimera', ['chrom1', 'strand1', 'start1', 'end1', 'qstart1', 'as1', 'chrom2', 'strand2', 'start2', 'end2', 'qstart2', 'as2', 'score', 'gap', 'ID', 'type', 'controltype'])

strand_conv = {True: '-', False: '+'}





##Functions

def get_best_mapping(choices, bestdistance):    
	bestchoice = max(choices, key =  lambda x: x.score)
	nonunique = filter(lambda x: (bestchoice.score-x.score)<bestdistance, choices)
	if(len(nonunique)==1):
		return bestchoice;
	else:
		return None;
	
	
def is_control(arw1, arw2):
	if(arw1.control or arw2.control):
		return 'control'
	else:
		return 'signal'
	
	
def _selected2chimera_forward(arw1, arw2, gap, type_):
	return Chimera(arw1.rname, strand_conv[arw1.aligned_read.is_reverse], arw1.aligned_read.reference_start, arw1.aligned_read.reference_end, arw1.qstart, arw1.AS, arw2.rname, strand_conv[arw2.aligned_read.is_reverse], arw2.aligned_read.reference_start, arw2.aligned_read.reference_end, arw2.qstart, arw2.AS, arw1.AS + arw2.AS, gap, arw1.qname, type_, is_control(arw1, arw2));

def _selected2chimera_reverse(arw1, arw2, gap, type_):
	return Chimera(arw2.rname, arw2.strand, arw2.aligned_read.reference_start, arw2.aligned_read.reference_end, arw1.qstart, arw2.AS, arw1.rname, arw1.strand, arw1.aligned_read.reference_start, arw1.aligned_read.reference_end, arw2.qstart, arw1.AS, arw1.AS + arw2.AS, gap, arw1.qname, type_, is_control(arw1, arw2));

def get_chimera(arw1, arw2, gap, funselect):
	if((arw1.rname == arw2.rname) and  (arw1.aligned_read.is_reverse == arw2.aligned_read.is_reverse)):
		ref_distance = arw2.aligned_read.reference_start - arw1.aligned_read.reference_start
		if(arw1.aligned_read.is_reverse):
			if(ref_distance > 0):
				return funselect(arw1, arw2, gap, 'csj')
			else:
				return funselect(arw1, arw2, gap, 'lsj')
		else:
			if(ref_distance > 0):
				return funselect(arw1, arw2, gap, 'lsj');
			else:
				return funselect(arw1, arw2, gap, 'csj');
	else:
		return funselect(arw1, arw2, gap, 'inter');

		
def arws2chimeras(arws, maxoverlap, maxgap, funselect):
	chimeras = [];
	for arw1, arw2 in combinations(sorted(arws, key = lambda x: x.qstart), 2):
		gap = arw2.qstart - arw1.qend;
		if(gap <= maxgap and gap >= maxoverlap):
			chimeras.append(get_chimera(arw1, arw2, gap, funselect))
	return chimeras;


def confirm_csj(csj, another):
	if(type(another) == Chimera):
		c = another.chrom1 == csj.chrom1 and another.chrom2 == csj.chrom1 and another.strand1 == csj.strand1 and another.strand2 == csj.strand2 and min(csj.start1, csj.start2)<=min(another.start1, another.start2) and max(csj.end1, csj.end2)>=max(another.end1, another.end2);
		#print 'chimera', c, another.type
		return c
	else:
		c = csj.chrom1==another.rname and csj.strand1==another.strand and min(csj.start1, csj.start2)<=another.aligned_read.reference_start and max(csj.end1, csj.end2)>=another.aligned_read.reference_end
		#print 'single', c
		return c
	
	
def confirm_lsj(lsj, another, forward, maxdistance):
	if(type(another) == Chimera):
		same_ref = another.chrom1 == lsj.chrom1 and another.chrom2 == lsj.chrom1 and another.strand1 == lsj.strand1 and another.strand2 == lsj.strand2
		distance = another.start2 - lsj.start2
	
		
	

def confirm_chimera(bestmates, maxdistance):
	if(type(bestmates[0]) == Chimera):
		if(bestmates[0].type == 'csj'):
			if(confirm_csj(bestmates[0], bestmates[1])):
				return bestmates[0];
	
	if(type(bestmates[1]) == Chimera):
		if(bestmates[1].type == 'csj'):
			if(confirm_csj(bestmates[1], bestmates[0])):
				return bestmates[1];
			
	#if(type(bestmates[0]) == Chimera):
		#if(bestmates[0].type == 'lsj'):
			#if(confirm_lsj(bestmates[0], bestmates[1], bestmates[0].strand1=='+', maxdistance)):
				#return bestmates[0];
			
	#if(type(bestmates[1]) == Chimera):
		#if(bestmates[1].type == 'lsj'):
			#if(confirm_lsj(bestmates[1], bestmates[0], bestmates[0].strand1=='-', maxdistance)):
				#return bestmates[1];
	

def demultiplex(paired_arws, maxdistance, maxoverlap, maxgap, concreward, bestdistance):
	if(all(paired_arws)):
		bestmates = [];
		for carws, funselect in zip(paired_arws, (_selected2chimera_forward, _selected2chimera_reverse)):
			choices =  carws + arws2chimeras(carws, maxoverlap, maxgap, funselect);
			bestmates.append(get_best_mapping(choices, bestdistance));
			
		if(all(bestmates)):
			return confirm_chimera(bestmates, maxdistance);
		else:
			return None

	else:
		return None




###########################################################################################
#function to write found single and chimeric hits to gff or doublegff files 


def parse2gff(best1):
	if(type(best1) == Chimera and best1.controltype == 'signal'):
		i1 = construct_gff_interval(best1.chrom1, best1.start1, best1.end1, best1.controltype, score=str(best1.score), strand=best1.strand1, source='un', frame='.', attrs=[("ID", best1.ID), ('qstart', best1.qstart1), ('AS', best1.as1), ('gap', best1.gap), ('mtype', best1.type)])
		
		i2 = construct_gff_interval(best1.chrom2, best1.start2, best1.end2, best1.controltype, score=str(best1.score), strand=best1.strand2, source='un', frame='.', attrs=[("ID", best1.ID), ('qstart', best1.qstart2), ('AS', best1.as2), ('gap', best1.gap), ('mtype', best1.type)])
		
		return "".join((str(i1), str(i2)));
	else:
		return None
		
			

			
		
############################################################################################
#run analysis on paired-end mapping results


# get score function
score2function = {'as': 'as_score', 'as_qstart': 'as_qstart_score', 'as_qstart_entropy': 'as_qstart_entropy_score', 'as_qstart_rstart': 'as_qstart_rstart_score', 'as_qstart_rstart_entropy': 'as_qstart_rstart_entropy_score'}
exec("from nrlbio.samlib import %s as key_score" % score2function[args.score]);

# get output file names 
fnames = [os.path.join(args.output, "%s.%s.gff" % (args.name, x)) for x in ['unique', 'control', 'unique_chimera', 'control_chimera']]

# fix overlap value 
fixed_maxoverlap = -args.maxoverlap


with open(fnames[0], 'w') as hsignal, open(fnames[1], 'w') as hcontrol, open(fnames[2], 'w') as hchimera_unique, open(fnames[3], 'w') as hchimera_control:
	types2handlers = { ('conc', 'signal'): hsignal, ('conc', 'control'): hcontrol, ('unconc', 'signal'): hchimera_unique, ('unconc', 'control'): hchimera_control }
	
	for rcount, paired_arws in enumerate(izip(generator_segments(args.path[0], key_score), generator_segments(args.path[1], key_score, secondmate=True))):
		best1 = demultiplex(paired_arws, args.maxdistance, fixed_maxoverlap, args.maxgap, args.concreward, args.bestdistance)
		ols = parse2gff(best1)
		if(ols):
			sys.stdout.write(ols) 
		
		
		
		#arws, score, mappingtype, swap, controltype = demultiplex(paired_arws, args.maxdistance, args.maxoverlap, args.maxgap, args.concreward, args.bestdistance)
		#if(arws):
			#types2handlers[(mappingtype, controltype)].write(parse2bed(arws, score, mappingtype, swap, controltype));
			
		if(rcount and rcount % 1000000 == 0):
			#sys.exit();
			sys.stderr.write("reads processed: %d\n" % rcount);
	
	

		
		
		
		
		