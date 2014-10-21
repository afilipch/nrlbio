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
#parser.add_argument('-n', '--norc', nargs = '?', required = True, type = str, help = "mapping was done not against genome");
#parser.add_argument('-f', '--filters', nargs = '+', required = True, type = str, help = "list of filters to apply");
args = parser.parse_args();

samfile = pysam.Samfile(args.path)

def _iteration(arlist):
	chimeras =  arlist2chimera(arlist, samfile, gap = 1, overlap = 4, score_function = chimera.as_gap_score)
	if(len(chimeras) == 1):
		print chimeras[0];	


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
			
			