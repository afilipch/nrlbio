#! /usr/lib/python
'''produce chimeras from merged and already filtered sam file'''
import argparse;
import os;

import pysam;
from nrlbio.filters_for_sam import *
from nrlbio.chimera import arlist2chimera
from nrlbio import chimera

parser = argparse.ArgumentParser(description='produce chimeras from merged and already filtered sam file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('--overlap', nargs = '?', default = 100, type = int, help = "max overlap allowed")
parser.add_argument('--gap', nargs = '?', default = 100, type = int, help = "max gap allowed")
args = parser.parse_args();

samfile = pysam.Samfile(args.path)

def _iteration(arlist):
	chimeras =  arlist2chimera(arlist, samfile, gap = args.gap, overlap = args.overlap, score_function = chimera.as_gap_score)
	if(len(chimeras) == 1):
		print chimeras[0].doublebed();	


arlist = [];
current_name = '';
for ar in samfile.fetch(until_eof=True):
	if(not ar.is_unmapped):
		
		if(current_name != ar.qname):
			_iteration(arlist);
				
			arlist = [ar];
			current_name = ar.qname;
			
		else:
			arlist.append(ar);
			
else:
	_iteration(arlist)
			
			