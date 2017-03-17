#! /usr/lib/python
'''Finds the positions of the binding sites on miRNA precursors'''

import argparse
import os
import sys
from collections import defaultdict

from pybedtools import BedTool, Interval

from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Finds the positions of the binding sites on miRNA precursors');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the miRNA targets, gff format ");
parser.add_argument("--precursors", nargs = '?', type=str, required=True, help="Path to the miRNA precursors, gff format")
args = parser.parse_args();


def get_postions(i1, i2):
	if(i1.strand == '-'):
		start = i1.end - i2.end;
	else:
		start = i2.start - i1.start;
	end = start+len(i2);
	return i2.attrs['Name'], start, end;


def list2interval(l):
	res = l[0]
	res.attrs['arm1'] = ":".join([str(x) for x in get_postions(res, l[1])]);
	if(len(l)>2):
		res.attrs['arm2'] = ":".join([str(x) for x in get_postions(res, l[2])]);
	return res;



###Get the intervals of the precursors augmented with the relative postions of a mature and a star miRNAs
precursors = []
curlist = [];
for interval in BedTool(args.precursors):
	if(interval[2] == "miRNA_primary_transcript"):
		if(curlist):
			precursors.append(list2interval(curlist))
		curlist = [interval];
	else:
		curlist.append(interval);
else:
	precursors.append(list2interval(curlist))
	

###Intersect miRNA tagerts recovered by ChiFlex with the precursor region and report miRNA;miRNA interactions;
def gff2dict(gff):
	return dict([tuple(x.split('=')) for x in gff.split(';')])

offset = 9;
for ii in BedTool(args.path).intersect(precursors, s = True, wo = True, f = 0.75):
	pattern, mirid, n_uniq = ii.attrs['pattern'], ii.attrs['mirid'], ii.attrs['n_uniq']
	pd = gff2dict(ii[2*offset-1])
	
	start_prec, end_prec = int(ii[offset+3]), int(ii[offset+4])
	start_target, end_target = int(ii[3]), int(ii[4])
	
	if(ii.strand == '-'):
		start = end_prec - end_target
	else:
		start = start_target - start_prec;
	end = start+end_target-start_target;
	
	target = "%s:%d:%d" % (mirid, start, end)
	
	res = construct_gff_interval(ii[offset], start_prec, end_prec, 'mirna_precursor', score='0', strand=ii.strand, source='un', frame='.', attrs=[('Name', pd['Name']), ('arm1', pd['arm1']), ('arm2', pd.get('arm2', '')), ('bsite', target), ('n_uniq', n_uniq), ('pattern', pattern)] )
	
	sys.stdout.write(str(res));
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
