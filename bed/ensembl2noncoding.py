#! /usr/bin/python
'''Extracts all noncoding regions from the ENSEMBL annotation'''


import sys;
import argparse;
from collections import defaultdict
import copy

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Extracts all noncoding regions from the ENSEMBL annotation');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the ensembl annotation file, gff3 format");
args = parser.parse_args();

transcripts = defaultdict(list)
for interval in BedTool(args.path):
	parent = interval.attrs.get('Parent')
	if(parent):
		ttype, name = parent.split(":")
		if(ttype == 'transcript'):
			transcripts[name].append(interval);
			

def transcript2noncoding(intervals):
	exons = [x for x in intervals if x[2] == 'exon']
	#exons.sort(key=lambda x: x.start)
	cds = [x for x in intervals if x[2] == 'CDS'][0]
	
	for exon in exons: 
		if(exon.start<cds.start and exon.end>cds.end):
			i1 = copy.copy(exon)
			i1.end = cds.start
			sys.stdout.write(str(i1))
			exon.start = cds.end
			sys.stdout.write(str(exon))
		elif(exon.start<cds.start):
			exon.end = min(cds.start, exon.end);
			sys.stdout.write(str(exon))
		elif(exon.end>cds.end):
			exon.start = max(exon.start, cds.end);
			sys.stdout.write(str(exon))
			
			
for name, intervals in transcripts.iteritems():
	types = set([x[2] for x in intervals])
	if('CDS' in types):
		transcript2noncoding(intervals)
	else:
		for interval in intervals:
			sys.stdout.write(str(interval))
	

