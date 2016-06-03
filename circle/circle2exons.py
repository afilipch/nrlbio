#! /usr/bin/python
'''Provides exonic structure of given circles on basis of transcript annotation(gff3 file)'''


import sys;
import argparse
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.pybedtools_extension import gff2bed


parser = argparse.ArgumentParser(description='Provides exonic structure of given circles on basis of transcript annotation(gff3 file)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the circles, bed/gff format");
parser.add_argument('--gff3', nargs = '?', required = True, type = str, help = "Path to the genome system annotation, gff3 format");
args = parser.parse_args();

def get_exons(interval, exons):
	'''Outputs all exons which may form a structure of the given interval(circle)
	
	interval pybedtools.Interval: genomic interval representing a circle (boundaries of the transcript)
	exons list: each element is a tuple representing genomic interval of an exon
	'''
	
	exons = sorted(list(set(exons)), key = lambda x: x[1])
	exons = filter(lambda x: int(x[1]) >= interval.start and int(x[2]) <= interval.end, exons)
	
	if(exons and interval.start == int(exons[0][1]) and  interval.end == int(exons[-1][2])):
		for exon in exons:
			exon = list(exon)
			exon[3] = interval.name
			exon[4] = '1'
			print "\t".join(exon);
	else:
		print "%s\t%d\t%d\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.end, interval.name, '0', interval.strand) 
		
	return None



#get exons from annotation gff3 file
exons = [];
for interval in BedTool(args.gff3):
	if('ID' in interval.attrs and interval.attrs['ID'].split(':')[0] == 'gene'):
		curname = interval.attrs['gene_id']
		enames = set()
	if(interval[2] == 'exon'):
		if(interval.name not in enames):
			enames.add(interval.name)
			interval.name = curname
			exons.append(gff2bed(interval))
		

	
#Get an intersection between circles and expms
bed = BedTool(args.path);
offset  = bed.field_count();
intersection = bed.intersect(b=exons, s=True, wao=True);

curname = ''
cexons = []
for interval in intersection:
	if(curname == interval.name):
		cexons.append(tuple(interval[offset:offset+6]))
	else:
		if(curname):
			get_exons(cinterval, cexons)
		cinterval = interval
		cexons = [tuple(interval[offset:offset+6])]
		curname = interval.name
else:
	get_exons(cinterval, cexons)
		

	


