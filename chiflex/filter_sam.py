#! /usr/lib/python
'''applies filter(s) to the given sam records'''
import argparse
import os;

import pysam;
from nrlbio.filters_for_sam import *

parser = argparse.ArgumentParser(description='applies filter(s) to the given sam records');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', required = True, type = str, help = "path to output filtered sam file");
parser.add_argument('-f', '--filters', nargs = '+', required = True, type = str, help = "list of filters to apply");
args = parser.parse_args();

samfile = pysam.Samfile(args.path)
filtered = pysam.Samfile(args.output, "wb", template=samfile)

for aligned_read in samfile.fetch(until_eof=True):
	if(not aligned_read.is_unmapped):
		for f in args.filters:
			if(not eval('%s(aligned_read)' % f)):
				break;
		else:
			filtered.write(aligned_read);
			
	