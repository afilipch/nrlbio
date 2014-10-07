#! /usr/lib/python
'''applies filter(s) to the given sam records'''
import argparse
import sys;

import pysam;
from nrlbio.filters_for_sam import *

parser = argparse.ArgumentParser(description='applies filter(s) to the given sam records');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', required = True, type = str, help = "path to output filtered sam file");
parser.add_argument('-f', '--filters', nargs = '+', required = True, choices=['repetitive', 'chimera_left', 'contain_repetitive'], type = str, help = "list of filters to apply");
args = parser.parse_args();

samfile = pysam.Samfile(args.path)
filtered = pysam.Samfile(args.output, "wb", template=samfile)

passed, removed = 0,0

for aligned_read in samfile.fetch(until_eof=True):
	if(not aligned_read.is_unmapped):
		for f in args.filters:
			if(not eval('%s(aligned_read)' % f)):
				removed +=1;
				break;
		else:
			filtered.write(aligned_read);
			passed +=1
			
			
sys.stderr.write("total reads: %d\npassed reads: %d\nremoved reads: %d\n" % (passed + removed, passed, removed));			
			
	