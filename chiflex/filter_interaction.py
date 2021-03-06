#! /usr/lib/python
'''applies filter(s) to the given interactions/chimeras'''
import argparse
import sys;

from nrlbio.generators import generator_doublebed
from nrlbio import filters_for_bed 
from nrlbio.formatting import string2fargs

parser = argparse.ArgumentParser(description='applies filter(s) to the given interactions/chimeras');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to interactions/chimeras file");
parser.add_argument('-f', '--filters', nargs = '+', required = True, choices=['distance', 'itype'], type = str, help = "list of filters to apply");
parser.add_argument('-a', '--arguments', nargs = '+', required = False, type = str, help = "key arguments for filters in format \'min_entropy=1.5,length=18 minright=15\'| None if no arguments provided for the filter. arguments have to be provided corresponding to the order of filters");
args = parser.parse_args();

_str2func={ 'distance': filters_for_bed.distance, 'itype': filters_for_bed.itype}
arguments = []
for a in args.arguments:
	if(a=='None'):
		arguments.append('');
	else:
		arguments.append(a)	

	
fargs = string2fargs(args.filters, arguments, _str2func)
passed, removed = 0,0

for doublebed in generator_doublebed(args.path):
	for fun, kwargs in fargs.items():
		if(not fun(doublebed,**kwargs)):
			removed +=1;
			break;		
	else:
		for interval in doublebed:
			sys.stdout.write(str(interval));
		passed +=1
			
			
sys.stderr.write("total interactions: %d\npassed interactions: %d\nremoved interactions: %d\n" % (passed + removed, passed, removed));	