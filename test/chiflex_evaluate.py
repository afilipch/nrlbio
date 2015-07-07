#! /usr/lib/python
'''Assignes to each read, if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits'''
import argparse
import os;
import sys;
from collections import defaultdict, namedtuple
from itertools import chain

import pysam;
from Bio import SeqIO;

from nrlbio.numerictools import overlap as get_overlap


parser = argparse.ArgumentParser(description='evaluates chiflex performance on artificially generated reads');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to initial reads. Fasta file");
parser.add_argument('-su', '--single_unique', nargs = '?', required = True, type = str, help = "Path to the unique single hits. Sam/bam format");
args = parser.parse_args();

Chimera = namedtuple('Chimera', 'typename chrom1 strand1 start1 stop1 ttype1 chrom2 strand2 start2 stop2 ttype2 left right gap')
Single = namedtuple('Single', 'typename chrom strand start stop ttype left right')
Shuffled = namedtuple('Shuffled', 'typename chrom strand start stop ttype length')
Random = namedtuple('Random', 'typename length');

#get initial reads
#_____________________________________________________________________________________________________________________________
def name2chimera(region, ttype, gaps):
	chrom1, strand1, start1, stop1, chrom2, strand2, start2, stop2= chain(*[x.split(":") for x in region.split("&")]);
	start1, stop1, start2, stop2 = [int(x) for x in (start1, stop1, start2, stop2)]
	
	ttype1, ttype2 = ttype.split(":");
	left, right, gap = [int(x) for x in gaps.split(":")]
	
	return Chimera('chimera', chrom1, strand1, start1, stop1, ttype1, chrom2, strand2, start2, stop2, ttype2, left, right, gap)
	
	
def name2single(region, ttype, gaps):
	chrom, strand, start, stop = region.split(":");
	start, stop = int(start), int(stop)

	left, right = [int(x) for x in gaps.split(":")]
	
	return Single('single', chrom, strand, start, stop, ttype, left, right)
	
	
def name2shuffled(region, ttype, gaps):
	chrom, strand, start, stop = region.split(":");
	start, stop = int(start), int(stop)

	length = int(gaps)
	
	return Shuffled('shuffled', chrom, strand, start, stop, ttype, length)	
	
	
def name2random(gaps):
	length = int(gaps)	
	return Random('random', length)		
	

reads = [];

def get_source(name):
	number, rtype, region, ttype, gaps = name.split("|")
	if(rtype == "chimera"):
		return name2chimera(region, ttype, gaps)
	elif(rtype == "single"):
		return name2single(region, ttype, gaps)
	elif(rtype == 'shuffled'):
		return name2shuffled(region, ttype, gaps)
	elif(rtype == 'random'):
		return name2random(gaps)
	else:
		sys.exit('Unrecognized read type. Has to be chimera|single|shuffled|random')
initial = defaultdict(int);

for i, seqrecord in enumerate(SeqIO.parse(args.path, 'fasta')):
	reads.append(get_source(seqrecord.name));
	
	
	
#get mapped reads
#_____________________________________________________________________________________________________________________________	
	
	
Comparison = namedtuple('Comparison', 'num length read_type mapped_type overlap');
strand_conv = {'-': True, '+': False}

def single_comparison(rnum, single, segment, segtype, segchrom):
	#for kv in vars(single).items():
		#print "%s\t%s" % tuple([str(x) for x in kv])
	if(single.typename == segtype):
		if(single.chrom == segchrom and strand_conv[single.strand] == segment.is_reverse):
			s, e = get_overlap((single.start, single.stop), (segment.reference_start, segment.reference_end));
			overlap = max(e-s, 0) 
			return Comparison(rnum, single.stop-single.start, single.typename, segtype, overlap);
		else:
			return Comparison(rnum, single.stop-single.start, single.typename, segtype, 0);
	else:
		return Comparison(rnum, 0, single.typename, segtype, 0);
		
		
def control_comparison(rnum, control, segment, segchrom):
	Comparison(rnum, control.length, single.typename, segtype, 0);
			
			
			
	
	
samfile = pysam.Samfile(args.single_unique);
for segment in samfile.fetch(until_eof=True):
	if(not segment.is_unmapped):
		rnum = int(segment.query_name.split("|")[0])
		segchrom = samfile.getrname(segment.tid)
		segtype = 'single'

		if(reads[rnum].typename == segtype):
			comparison = single_comparison(rnum, reads[rnum], segment, segtype, segchrom)
			print segment.query_name, segchrom, segment.is_reverse, segment.reference_start, segment.reference_end
			print
			print reads[rnum]
			print
			print comparison
			print "_"*120
	
	
	