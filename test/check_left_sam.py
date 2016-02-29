#! /usr/lib/python
'''Checks sensetivity and specificity of the left(miRNA) chimeric part mapping/demultiplexing'''
import argparse
import os;
import sys
from collections import defaultdict, namedtuple
from itertools import product, combinations, izip

import pysam;
from pybedtools import BedTool


from nrlbio.samlib import ArWrapper, demultiplex_read_hits
from nrlbio.generators import generator_segments




parser = argparse.ArgumentParser(description='Checks sensetivity and specificity of the left(miRNA) chimeric part mapping/demultiplexing');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the doublebed file of artificially generated chimeras");
parser.add_argument('--sam', nargs = '?', required = True, type = str, help = "Path to the mapped reads sam/bam format");
parser.add_argument('--overlap', nargs = '?', default = 12, type = int, help = "Min overlap required for mapping to be counted as correct");
parser.add_argument('--nonunique', nargs = '?', default = None, type = str, help = "Path to the auxillary file of nonunique region. If set hit will be counted as true positive if it source is one of the identical regions");
args = parser.parse_args();

nu_dict = defaultdict(list);
if(args.nonunique):
	for interval in BedTool(args.nonunique):
		nu_dict[int(interval.name)].append(interval);
		


generated = [];
last_num = ''

for interval in BedTool(args.path):
	if(interval.name != last_num):
		generated.append(interval);
	last_num = interval.name;


def arw_vs_interval(arw, interval, overlap):
	if(arw.rname == interval.chrom and arw.strand == interval.strand):
		return min(interval.end, arw.aligned_read.reference_end) - max(interval.start, arw.aligned_read.reference_start) >= overlap;
	else:
		return False
	
	
def interval_vs_nonunique(interval, nonunique, overlap):
	for nu in nonunique:
		if(min(interval.end, nu.end) - max(interval.start, nu.start) >= overlap):
			return True
	else:
		return False;
		

def simple_check(arws, interval, overlap, nu_dict):
	return any([arw_vs_interval(arw, interval, overlap) for arw in arws])

def nonunique_check(arws, interval, overlap, nu_dict):
	if(arw_vs_interval(arws[0], interval, overlap)):
		return True
	else:
		return interval_vs_nonunique(interval, nu_dict[arws[0].aligned_read.get_tag('CI')], overlap)

if(args.nonunique):
	check = nonunique_check
else:
	check = simple_check


checked = {};
for arws in generator_segments(args.sam):
	idnum = int(arws[0].qname.split("|")[0])
	checked[idnum] = int(check(arws, generated[idnum], args.overlap, nu_dict));
	
	
mirtypes = set(['mirna_chimera', 'mirna_single', 'mirna_shuffled', 'mirna_random'])
mtypes = {1: 'mapped', 0: 'fmapped', -1: 'unmapped'}

for idnum, interval in enumerate(generated):
	if(interval.attrs['rtype'] in mirtypes):
		interval.attrs['mapped'] = mtypes[checked.get(idnum, -1)]
		sys.stdout.write(str(interval))
	
	
	