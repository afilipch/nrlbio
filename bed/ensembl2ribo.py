#! /usr/bin/python
'''Extracts all ribosomal regions from the ENSEMBL annotation'''


import sys;
import argparse;
from collections import defaultdict
import copy

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Extracts all ribosomal regions from the ENSEMBL annotation');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the ensembl annotation file, gff3 format");
args = parser.parse_args();

for interval in BedTool(args.path):
	btype = interval.attrs.get('biotype', 'none')
	ttype = interval.attrs.get('ID', 'none').split(':')[0]
	if(ttype=='gene' and btype in ('rRNA', 'Mt_rRNA')):
		sys.stdout.write(str(interval))