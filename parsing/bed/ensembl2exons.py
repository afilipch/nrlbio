#! /usr/bin/python
'''Extracts exons from a given ensembl annotation gff3 file'''


import sys;
import argparse

from pybedtools import BedTool

from nrlbio.pybedtools_extension import gff2bed


parser = argparse.ArgumentParser(description='Extracts exons from a givne ensembl annotation gff3 file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome system annotation, gff3 ensembl format");
args = parser.parse_args();


exons = [];
for interval in BedTool(args.path):
	if('ID' in interval.attrs and interval.attrs['ID'].split(':')[0] == 'gene'):
		curname = interval.attrs['gene_id']
		enames = set()
	if(interval[2] == 'exon'):
		if(interval.name not in enames):
			enames.add(interval.name)
			interval.name = curname
			sys.stdout.write(str(gff2bed(interval)))