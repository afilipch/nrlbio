#! /usr/bin/python
'''Provides exonic structure of given circles on basis of transcript annotation(gff3 file) in bed12 format'''


import sys;
import argparse
from collections import defaultdict

from pybedtools import BedTool
import pysam


parser = argparse.ArgumentParser(description='Extracts read coverage for the given intervals');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the interval, bed/gff format");
parser.add_argument('--sortedbam', nargs = '?', required = True, type = str, help = "Path to the sorted read hits, bam format");
parser.add_argument('--collapsed', nargs = '?', default = False,  const = True, type = bool, help = "If set, reads are assumed to be collapsed with collapse.pl script");
parser.add_argument('--totalcoverage', nargs = '?', type = int, help = "If set, total covergae will be not recalculated");
args = parser.parse_args();

stranddict = {'-': True, '+': False}


def segment2count_noncollapsed(segment):
	return segment.reference_length

def segment2count_collapsed(segment):
	return segment.reference_length*int(segment.query_name.split('_')[-1])


def column2count_noncolllapsed(pileupcolumn, strand):
	c = 0
	for read in pileupcolumn.pileups:
		if(read.alignment.is_reverse == strand):
			c += 1
	return c

def column2count_collapsed(pileupcolumn, strand):
	c = 0
	for read in pileupcolumn.pileups:
		if(read.alignment.is_reverse == strand):
			c += int(read.alignment.query_name.split('_')[-1])
	return c

if(args.collapsed):
	segment2count = segment2count_collapsed
	column2count = column2count_collapsed
else:
	segment2count = segment2count_noncollapsed
	column2count = column2count_noncolllapsed
	
	

samfile = pysam.AlignmentFile(args.sortedbam, 'rb')

intervals = [];
for n, interval in enumerate(BedTool(args.path)):
	count = 0;
	strand = stranddict[interval.strand];
	for pileupcolumn in samfile.pileup(interval.chrom, interval.start, interval.end):
		count += column2count(pileupcolumn, strand)
	interval.attrs['raw_coverage'] = str(count)
	intervals.append(interval)
	if(n and n % 1000 == 0):
		sys.stderr.write("%d intervals are processed\n" % n)
	
samfile.close();


if(args.totalcoverage):
	totalcoverage = args.totalcoverage
else:
	samfile = pysam.AlignmentFile(args.sortedbam, 'rb');
	totalcoverage = 0;
	for segment in samfile.fetch(until_eof=True):
		totalcoverage += segment2count(segment);
	samfile.close();
	
sys.stderr.write("%d is a total coverage\n" % totalcoverage)



for interval in intervals:
	interval.attrs['norm_coverage'] = "%1.7f" % ((float(interval.attrs['raw_coverage'])*1000000)/(totalcoverage*len(interval)))
	sys.stdout.write(str(interval))

	
#for read in samfile.fetch('chr2', 19500000, 19505000):
	#print "\t".join(str(read).split("\t")[:5])
	
	
#print 
#print

#for read in samfile.fetch('chr2', 19500000, 19505000):
	#print "\t".join(str(read).split("\t")[:5])
	
#for pileupcolumn in samfile.pileup('chr2', 19500000, 19505000):
	#print pileupcolumn.pos
	#for read in pileupcolumn.pileups:
		#print read.alignment.query_name
	#print
	
	
