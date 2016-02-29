#! /usr/lib/python
'''creates bed12 file from bed file with exons. Exons of the same transcript should have the same id''' 
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool, Interval;

parser = argparse.ArgumentParser(description='creates bed12 file from bed file with exons. Exons of the same transcript should have the same id');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "bed file with exons. Exons within the same transcript must have the same id");
#parser.add_argument('-f', '--fasta', nargs = '?', required = True, type = str, help = "path to fasta file(genome)");
args = parser.parse_args();



def fetch_bed12(name, intervals):
	intervals.sort(key = lambda x: x.start);
	chrom = intervals[0].chrom
	start = str(intervals[0].start)
	end = str(intervals[-1].stop)
	score = '0'
	strand = intervals[0].strand
	
	thick_start = start
	thick_end = end
	itemRGB = '0,0,255'
	
	block_count = str(len(intervals))
	block_sizes = ",".join([str(len(i)) for i in intervals])
	block_starts = ",".join([str(i.start) for i in intervals])
	print "\t".join([chrom, start, end, name, score, strand, thick_start, thick_end, itemRGB, block_count, block_sizes, block_starts])
	

curchrom = '';
chrdict = defaultdict(list);
for i in BedTool(args.path).sort():
	if(chrdict and i.chrom != curchrom):
		for k, v in chrdict.iteritems():
			fetch_bed12(k, v);
		chrdict = defaultdict(list);
		chrdict[i.name] = [i];
		curchrom = i.chrom
	else:
		chrdict[i.name].append(i);
else:
	for k, v in chrdict.iteritems():
		fetch_bed12(k, v);