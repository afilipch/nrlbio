#! /usr/lib/python
'''filters gff intervals on basis of their attributes values'''
import argparse;
import sys;
import operator

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='filters gff intervals on basis of their attributes values');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to interactions/chimeras file");
parser.add_argument('-a', '--attributes', nargs = '+', default=[], type = str, help = "list of the attributes to filter with");
parser.add_argument('-c', '--cutoffs', nargs = '+', default=[], type = str, help = "list of the cutoffs to be applied for each attribute correspondingly. The order of cutoffs has to be consistent with \'--attributes\'");
parser.add_argument('-m', '--modes', nargs = '+', default=[], choices=['less', 'equalless', 'greater', 'equalgreater', 'equal', 'nonequal'], type = str, help = "list of the filtering modes to be applied for each attribute and cutoff correspondingly. The order of cutoffs has to be consistent with \'--attributes\'");
args = parser.parse_args();


if(len(args.attributes) == len(args.cutoffs) == len(args.modes)):
	pass;
else:
	raise ValueError("Numbers of arguments provided to the \'--attributes\', \'--cutoffs\', \'--modes\' have to be equal")

	
mode2fun = {"less": operator.lt, "equalless": operator.le, 'greater': operator.gt, "equalgreater": operator.ge, 'equal': operator.eq, "nonequal": operator.ne}	

def _filter(interval, rules):
	for attr, fun, cutoff, type_ in rules:
		if(fun(type_(interval.attrs[attr]), cutoff)):
			pass;
		else:
			return False;
	else:
		return True


#get rules
rules = [];
for attr, mode, cutoff in zip(args.attributes, args.modes, args.cutoffs):
	if(cutoff.lstrip('-').isdigit()):
		c = int(cutoff)
		type_ = int
	else:
		try:
			c = float(cutoff);
			type_ = float
		except:
			c = str(cutoff);
			type_ = str
			
	rules.append((attr, mode2fun[mode], c, type_));
	
	
passed, removed = 0,0
for interval in BedTool(args.path):
	if(_filter(interval, rules)):
		sys.stdout.write(str(interval))
		passed+=1;
	else:
		removed+=1;
			
			
sys.stderr.write("total intervals: %d\npassed intervals: %d\nremoved intervals: %d\n" % (passed + removed, passed, removed));		
