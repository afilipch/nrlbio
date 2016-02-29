#! /usr/lib/python
'''Demultiplexing of raw mapping reasults for paired-end sequncing'''
import argparse
import os;
import sys
from collections import defaultdict
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
parser.add_argument('--bestdistance', nargs = '?', default = 10, type = float, help = "minimal distance allowed between the best and the second best hit. If the actual distance is less, than hit will be assigned as nonunique");
parser.add_argument('--maxdistance', nargs = '?', default = 1000, type = int, help = "Max distance allowed for paired-end hits to be considered as concordant alignment");
parser.add_argument('--maxgap', nargs = '?', default = 0, type = int, help = "Max gap between hits on the same read allowed to consider this read(one mate) to be chimeric");
parser.add_argument('--maxoverlap', nargs = '?', default = 6, type = int, help = "Max overlap between hits on the same read allowed to consider this read(one mate) to be chimeric");
parser.add_argument('--concreward', nargs = '?', default = 2, type = float, help = "Multiplier for the score for concordant mappings");
parser.add_argument('--verbose', nargs = '?', default = 0, const = 1000000, type = int, help = "If set, progress of the demultiplexing will be reported");
args = parser.parse_args();

# parse [score] argument
score2function = {'as': 'as_score', 'as_qstart': 'as_qstart_score', 'as_qstart_entropy': 'as_qstart_entropy_score', 'as_qstart_rstart': 'as_qstart_rstart_score', 'as_qstart_rstart_entropy': 'as_qstart_rstart_entropy_score'}
exec("from nrlbio.samlib import %s as key_score" % score2function[args.score]);

# get output file names 
fnames = [os.path.join(args.output, "%s.%s.gff" % (args.name, x)) for x in ['unique', 'control', 'unique_chimera', 'control_chimera']]


def clear_dict(d):
	for v in d.values():
		v[:]=[];
		

def doublet_score(arws):
	return sum([x.score for x in arws])

def triplet_score(ch1, ch2, arw):
	return ch1.score + ch2.score + arw.score;
		
		
def remove_nonunique(mates, bestdistance):    
	bestmate = max(mates, key =  lambda x: x[1])
	#return bestmate
	nonunique = filter(lambda x: (bestmate[1]-x[1])<bestdistance, mates)
	if(len(nonunique)==1):
		return bestmate;
	else:
		return None;
	
def is_control(arws):
	if(any([x.control for x in arws])):
		return 'control'
	else:
		return 'signal'


		
def mate2chimera(arws, maxoverlap, maxgap):
	distance = max([x.qstart for x in arws]) - min([x.qend for x in arws])
	if(-1*maxoverlap<=distance<=maxgap):
		return tuple(sorted(arws, key = lambda x: x.qstart))
	else:
		return None
	
def triplet2chimera(matechimeras, arws, maxdistance, concreward, swap):
	triplets = [];
	for (ch1,ch2), arw in product(matechimeras, arws):
		if(ch2.rname == arw.rname):
			distance = max(ch2.aligned_read.reference_start, arw.aligned_read.reference_start) - min(ch2.aligned_read.reference_end, arw.aligned_read.reference_end)
			if(distance<=maxdistance):
				triplets.append(((ch1, ch2, arw), triplet_score(ch1, ch2, arw)*concreward, 'unconc', swap));
	return triplets;
	
	
def doublet2chimera(arws, maxdistance, concreward):
	if(arws[0].rname == arws[1].rname):
		distance = max([x.aligned_read.reference_start for x in arws]) - min([x.aligned_read.reference_end for x in arws])
		if(distance<=maxdistance):
			return (arws, doublet_score(arws)*concreward, 'conc', False);
		else:
			return (arws, doublet_score(arws), 'unconc', False);
	else:
		return (arws, doublet_score(arws), 'unconc', False);
	

def demultiplex(paired_segments, maxdistance, maxoverlap, maxgap, concreward, bestdistance):
	mates = [];
	for arws in product(paired_segments[0], paired_segments[1]):
		mates.append(doublet2chimera(arws, maxdistance, concreward))
			
	matechimeras = defaultdict(list)
	for mate_num, all_arws in enumerate(paired_segments):
		for arws in combinations(all_arws, 2):
			mc = mate2chimera(arws, maxoverlap, maxgap)
			if(mc):
				matechimeras[mate_num].append(mc);
		
	mates += triplet2chimera(matechimeras[0], paired_segments[1], maxdistance, concreward, False)
	mates += triplet2chimera(matechimeras[1], paired_segments[0], maxdistance, concreward, True)
	
	if(mates):
		bestmate = remove_nonunique(mates, bestdistance);
		if(bestmate):
			return bestmate[0], bestmate[1], bestmate[2], bestmate[3], is_control(bestmate[0])
			
	return [None]*5




###########################################################################################
#function to write found single and chimeric hits to gff or doublegff files 

strand_conv = {True: '-', False: '+'}

def parse2bed(arws, score, mappingtype, swap, controltype):
	if(mappingtype == 'conc'):
		arw = arws[0]
		start = min([x.aligned_read.reference_start for x in arws])
		stop = max([x.aligned_read.reference_end for x in arws])
		strand = strand_conv[arw.aligned_read.is_reverse]
		AS = sum([x.AS for x in arws])
		qstart = min([x.qstart for x in arws])
		return str(construct_gff_interval(arw.rname, start, stop, controltype, score=str(score), strand=strand, source='un', frame='.', attrs=[("ID", arw.qname.split(":")[0]), ('qstart', qstart), ('AS', AS)]))

	else:
		if(len(arws)==2):
			intervals = []
			strand = strand_conv[arws[0].aligned_read.is_reverse];
			for arw in arws:
				intervals.append(construct_gff_interval(arw.rname, arw.aligned_read.reference_start, arw.aligned_read.reference_end, controltype, score=str(score), strand=strand, source='un', frame='.', attrs=[("ID", arw.qname), ('qstart', arw.qstart), ('AS', arw.AS), ('gap', 0), ('matech', 'doublet')]))
			return "".join([str(x) for x in intervals])	

		else:
			gap = max([x.qstart for x in arws[:2]]) - min([x.qend for x in arws[:2]])
			if(swap):
				strand1 = strand_conv[not arws[0].aligned_read.is_reverse]
				strand2 = strand_conv[arws[2].aligned_read.is_reverse]
			else:
				strand1 = strand_conv[arws[0].aligned_read.is_reverse]
				strand2 = strand_conv[not arws[2].aligned_read.is_reverse]
			
			arw = arws[0];
			i1 = construct_gff_interval(arw.rname, arw.aligned_read.reference_start, arw.aligned_read.reference_end, controltype, score=str(score), strand=strand1, source='un', frame='.', attrs=[("ID", arw.qname), ('qstart', arw.qstart), ('AS', arw.AS), ('gap', gap), ('matech', 'triplet')])
			#matenum =  arw.qname.split(":")[1]
			
			chrom = arws[2].rname
			start = min([x.aligned_read.reference_start for x in arws[1:]])
			stop = max([x.aligned_read.reference_end for x in arws[1:]])
			strand = strand_conv[arws[2].aligned_read.is_reverse]
			name = arws[2].qname
			AS = sum([x.AS for x in arws[1:]])
			qstart = arws[1].qstart
			
			i2 = construct_gff_interval(chrom, start, stop, controltype, score=str(score), strand=strand2, source='un', frame='.', attrs=[("ID", name), ('qstart', qstart), ('AS', AS), ('gap', gap), ('matech', 'triplet')])
			
			if(swap):
				i1, i2 = i2, i1
			
			return "".join((str(i1), str(i2)))
			
		
		

	#l = list([x.rname.split("|")[-1] for  x in arws])
	#return "\t".join(l)
	
#######################################test functions
	
def just_pairs(paired_segments):
	'''function for first step comparison'''
	#print '*'*140
	for p in product(*paired_segments.values()):
		l = [x.rname.split("|")[-1] for  x in p]
		if(l and l[0]!=l[1]):
			l.sort();
			print "\t".join(l);
	
	
def check_reverse(paired_segments, maxdistance, maxoverlap, maxgap, concreward):
	for arw in paired_segments['1']:
		print arw.aligned_read.is_reverse
#####################################################		
	
demultiplex = demultiplex
	


#samfile = pysam.Samfile(args.path);
#paired_segments = defaultdict(list);
#current_name = '';
#rcount = 0;
#read2mate = {True: '1', False: '2'}

with open(fnames[0], 'w') as hsignal, open(fnames[1], 'w') as hcontrol, open(fnames[2], 'w') as hchimera_unique, open(fnames[3], 'w') as hchimera_control:
	types2handlers = { ('conc', 'signal'): hsignal, ('conc', 'control'): hcontrol, ('unconc', 'signal'): hchimera_unique, ('unconc', 'control'): hchimera_control }
	
	for rcount, paired_segments in enumerate(izip(generator_segments(args.path[0], key_score), generator_segments(args.path[1], key_score))):
		arws, score, mappingtype, swap, controltype = demultiplex(paired_segments, args.maxdistance, args.maxoverlap, args.maxgap, args.concreward, args.bestdistance);
		if(arws):
			types2handlers[(mappingtype, controltype)].write(parse2bed(arws, score, mappingtype, swap, controltype));
			
		if(rcount and rcount % 1000000 == 0):
			#sys.exit();
			sys.stderr.write("reads processed: %d\n" % rcount);
	
	
	#for segment in samfile.fetch(until_eof=True):
		#if(not segment.is_unmapped):
			#rname = samfile.getrname(segment.tid)
			#arw = ArWrapper(segment, rname, score_function=key_score, add_nr_tag=False)
			#basename = segment.query_name.split(':')[0];
			#mate = read2mate[segment.is_read1];
			
			#if(current_name != basename):
				#rcount+=1
				#arws, score, mappingtype, controltype = demultiplex(paired_segments, args.maxdistance, args.maxoverlap, args.maxgap, args.concreward, args.bestdistance);
				#if(arws):
					#types2handlers[(mappingtype, controltype)].write(parse2bed(arws, score, mappingtype, controltype));
				#clear_dict(paired_segments);
				#current_name = basename;
				
			#paired_segments[mate].append(arw);
			
		#if(args.verbose and (not rcount % args.verbose)):
			#print "%d sam entries have been processed" % rcount
			#break;
			
	#else:
		#arws, score, mappingtype, controltype = demultiplex(paired_segments, args.maxdistance, args.maxoverlap, args.maxgap, args.concreward, args.bestdistance);
		#if(arws):
			#types2handlers[(mappingtype, controltype)].write(parse2bed(arws, score, mappingtype, controltype));
	
		
		
		
		
		