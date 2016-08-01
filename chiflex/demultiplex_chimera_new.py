#! /usr/lib/python
'''Assignes to each read, if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits'''
import argparse
import os;
import sys;

import pysam;

import nrlbio.statistics.sam as sam_statistics
from nrlbio.samlib import ArWrapper, demultiplex_read_hits
from nrlbio.chimera import arwlist2chimeras
#from nrlbio.chimera import demultiplex as demultiplex_ch
from nrlbio.generators import generator_segments



parser = argparse.ArgumentParser(description='Assignes to each read if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', default = "sam", type = str, help = "path to the output folder");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for output files, should reflect nature of mapping reference");

parser.add_argument('-mg', '--maxgap', nargs = '?', default = 4, type = int, help = "Max gap allowed between parts of chimera. It is also used to calculate chimera score. The more the maxgap, the less important is the gap.");
parser.add_argument('-mo', '--maxoverlap', nargs = '?', default = 4, type = int, help = "Max overlap allowed between parts of chimera. It is also used to calculate chimera score. The more the maxoverlap, the less important is the overlap.");
parser.add_argument('-sd', '--splice_distance', nargs = '?', default = 1000000, type = int, help = "If the distance between chimeric parts on a reference is less than a splice_distance, then this chimera will be preferentially selected as potential splice junction");

parser.add_argument('--s_distance', nargs = '?', default = 12, type = float, help = "minimal distance allowed between the best and the second best hit. If the actual distance is less, than hit will be assigned as nonunique");
parser.add_argument('--ch_distance', nargs = '?', default = 12, type = float, help = "minimal distance allowed between the best and the second best chimera. If the actual distance is less, than chimera will be assigned as nonunique");
args = parser.parse_args();


def select_splicing(chimeras, splice_distance):
	return [x for x in chimeras if (x.arws[0].rname == x.arws[1].rname and x.arws[0].strand == x.arws[1].strand and abs(x.arws[0].aligned_read.reference_start - x.arws[1].aligned_read.reference_start)<=splice_distance) ]

def chimera2score(chimera, gapbase, overlapbase):
	'''Score function for Chimera based on alignment score and gap'''
	if(chimera.gap>0):
		gap = chimera.gap
		overlap = 0
	else:
		gap = 0;
		overlap = -chimera.gap
	
	return chimera.AS*(1 - gap/gapbase)*(1-overlap/overlapbase)


def get_chimeras(arwlist, maxgap, maxoverlap, splice_distance):
	chimeras = arwlist2chimeras(arwlist, maxgap, maxoverlap)
	spliced = select_splicing(chimeras, splice_distance)
	if(spliced):
		chimeras = spliced;





#mappings derived from the same read are pulled together. Collapsed into one (or more, in case of non-unique mappings) wrapping read object. 
#then they are written to a new destination, according to their source: real, or control
for arwlist in generator_segments(args.path, key_score):
	get_chimeras(arwlist, args.maxgap, args.maxoverlap, args.splice_distance)
	#_iteration(arwlist, args.s_distance, args.ch_distance)

	

			

			
			