#! /usr/lib/python
'''applies filter(s) to the given interactions/chimeras'''
import argparse
import sys;

from nrlbio.generators import generator_doublebed
from nrlbio.filters_for_bed import *

parser = argparse.ArgumentParser(description='applies filter(s) to the given interactions/chimeras');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to interactions/chimeras file");
parser.add_argument('-f', '--filters', nargs = '+', required = True, choices=['distance', 'itype'], type = str, help = "list of filters to apply");
parser.add_argument('-a', '--arguments', nargs = '+', required = False, type = str, help = "key arguments for filters in format \'min_entropy=1.5,length=18 minright=15\'| None if no arguments provided for the filter. arguments have to be provided corresponding to the order of filters");
args = parser.parse_args();

#some parsing of agrs.arguments to make them ready to paste into eval statement
arguments = [];
if(args.arguments):
	if(len(args.arguments) == len(args.filters)):
		for a in args.arguments:
			if(a=='None'):
				arguments.append('');
			else:
				arguments.append("".join((",", a)))
	else:
		raise AttributeError('number of arguments should be equal to number of filters, or not provided at all')
else:
	arguments = ['']*len(args.filters);
	




passed, removed = 0,0

for doublebed in generator_doublebed(args.path):
	for f, a in zip(args.filters, arguments):
		if(not eval('%s(doublebed%s)' % (f, a))):
			removed +=1;
			break;
	else:
		for interval in doublebed:
			sys.stdout.write(str(interval));
		passed +=1
			
			
sys.stderr.write("total interactions: %d\npassed interactions: %d\nremoved interactions: %d\n" % (passed + removed, passed, removed));	