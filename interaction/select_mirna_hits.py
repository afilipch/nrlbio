#! /usr/bin/python
'''Selects mapping hits of chimeric target parts only for mirna ids provided''' 
import argparse
import sys;

import pysam;

from nrlbio.generators import generator_doublebed



parser = argparse.ArgumentParser(description='Selects mapping hits of chimeric target parts only for mirna ids provided');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the mirna:target interactions to select from, double gff");
parser.add_argument('--mirids', nargs = '+', required=True, type = str, help = "miRNA ids to select");
parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Table (rid2iid.tsv) which connects interaction ids with read ids");
parser.add_argument('--sam', nargs = '?', required=True, type = str, help = "Path to mapping hits, sam/bam file");
parser.add_argument('--output', nargs = '?', required = True, type = str, help = "Path for the selected sam/bam file");
parser.add_argument('--inverse', nargs = '?', default=False, const=True, type = bool, help = "Inverse selection. Selects all mapping hits of chimeric target parts, but for mirna ids provided");
args = parser.parse_args();

mirids = set(args.mirids)
iids = set()

if(args.inverse):
	for i1, i2 in generator_doublebed(args.path):
		if(i1.chrom not in mirids):
			iids.add(i1.name.split("|")[0])
else:		
	for i1, i2 in generator_doublebed(args.path):
		if(i1.chrom in mirids):
			iids.add(i1.name.split("|")[0]);
			
	
rids = set();

with open(args.table) as f:
	for l in f:
		rid, iid = l.strip().split("\t");
		if(iid in iids):
			rids.add(rid);
			
samfile = pysam.Samfile(args.sam);
selected = pysam.Samfile(args.output, "wb", template=samfile)

for segment in samfile.fetch(until_eof=True):
	if(segment.query_name in rids):
		selected.write(segment);
		
		
samfile.close();
selected.close();
