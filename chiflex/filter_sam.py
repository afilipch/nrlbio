#! /usr/bin/python
'''Applies filter(s) to the given sam records'''
import argparse
import sys;

import pysam;
from nrlbio.filters_for_sam import *
from nrlbio import filters_for_sam
from nrlbio.formatting import string2fargs

parser = argparse.ArgumentParser(description='Applies filter(s) to the given sam records');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', required = True, type = str, help = "Path to output filtered sam file");
parser.add_argument('-f', '--filters', nargs = '+', required = True, choices=['repetitive', 'chimera_right', 'contain_repetitive'], type = str, help = "List of filters to apply");
parser.add_argument('-a', '--arguments', nargs = '+', required = False, type = str, help = "Key arguments for filters in format \'min_entropy=1.5,length=18 minright=15\'| None if no arguments provided for the filter. arguments have to be provided corresponding to the order of filters");
args = parser.parse_args();



_str2func={ 'repetitive': filters_for_sam.repetitive, 'chimera_right': filters_for_sam.chimera_right, 'contain_repetitive': filters_for_sam.contain_repetitive}
arguments = []
if(args.arguments):
	for a in args.arguments:
		if(a=='None'):
			arguments.append('');
		else:
			arguments.append(a)	
else:
	arguments = ['']*len(args.filters)

	
fargs = string2fargs(args.filters, arguments, _str2func)
samfile = pysam.Samfile(args.path)
filtered = pysam.Samfile(args.output, "wb", template=samfile)

passed, removed = 0,0

for aligned_read in samfile.fetch(until_eof=True):
	if(not aligned_read.is_unmapped):
		for fun, kwargs in fargs.items():
			if(not fun(aligned_read,**kwargs)):
				removed +=1;
				break;		
		else:
			passed+=1
			filtered.write(aligned_read);

				
			
sys.stderr.write("total hits: %d\npassed hits: %d\nremoved hits: %d\n" % (passed + removed, passed, removed));	

	
			
	