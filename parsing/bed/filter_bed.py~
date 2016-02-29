'''filters bed/gff file according to some filters'''

import argparse
import sys;
from collections import defaultdict, Counter;

from pybedtools import BedTool;




parser = argparse.ArgumentParser(description='filters bed/gff file according to some filters');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed/gff file");
parser.add_argument('-f', '--filters', nargs = '+', required = True, choices=['length'], type = str, help = "list of filters to apply");
parser.add_argument('-a', '--arguments', nargs = '+', required = False, type = str, help = "key arguments for filters in format \'min_entropy=1.5,length=18 minright=15\'| None if no arguments provided for the filter. arguments have to be provided corresponding to the order of filters");
args = parser.parse_args();

#list of filters

def length(interval, minlength=None, maxlength=None):
	if(minlength and len(interval) < minlength):
		return False;
	elif(maxlength and len(interval) >= maxlength):
		return False;
	else:
		return True;


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
	
#print arguments	


passed, removed = 0,0

for interval in BedTool(args.path):
	for f, a in zip(args.filters, arguments):
		if(not eval('%s(interval%s)' % (f, a))):
			removed +=1;
			break;
	else:
		sys.stdout.write(str(interval));
		passed +=1
		
sys.stderr.write("total entries: %d\npassed entries: %d\nremoved entries: %d\n" % (passed + removed, passed, removed));			





