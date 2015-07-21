#! /usr/lib/python
'''produce chimeras from merged and already filtered sam file'''
import argparse;
import os;

import pysam;

from nrlbio.filters_for_sam import *
from nrlbio.chimera import Chimera, as_score
from nrlbio.samlib import ArWrapper
from nrlbio.generators import generator_doublesam

parser = argparse.ArgumentParser(description='produce chimeras from merged and already filtered sam file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('--overlap', nargs = '?', default = 100, type = int, help = "max overlap allowed")
parser.add_argument('--gap', nargs = '?', default = 100, type = int, help = "max gap allowed")
args = parser.parse_args();

samfile = pysam.Samfile(args.path)
for segments in generator_doublesam(samfile):
	segchroms = [samfile.getrname(segment.tid) for segment in segments]
	ar_wrappers = [ArWrapper(segment, segchrom) for segment, segchrom in zip(segments, segchroms)]
	ch = Chimera(ar_wrappers, as_score)
	print ch.doublebed()

#def _iteration(arwlist):
	#chimeras =  arwlist2chimera(arwlist, gap = args.gap, overlap = args.overlap, score_function = as_gap_score)
	#if(len(chimeras) == 1):
		#print chimeras[0].doublebed();	


#arwlist = [];
#current_name = '';
#for ar in samfile.fetch(until_eof=True):
	#if(not ar.is_unmapped):
		#rname = samfile.getrname(ar.tid)
		#arw = ArWrapper(ar, rname, add_nr_tag=False)
		
		#if(current_name != arw.qname):
			#_iteration(arwlist);
				
			#arlist[:] = [arw];
			#current_name = arw.qname;
			
		#else:
			#arlist.append(arw);
			
#else:
	#_iteration(arwlist)
			
			