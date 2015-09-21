#! /usr/lib/python
'''applies filter(s) to the given interactions/chimeras'''
import argparse
import sys;
import os;

from pybedtools import BedTool
from collections import defaultdict

parser = argparse.ArgumentParser(description='stratifies given gff file on basis of some attributes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to gff file");
parser.add_argument('-a', '--attributes', nargs = '+', required = True, type = str, help = "list of attributes used to stratify given gff entries");
parser.add_argument('-o', '--output', nargs = '?', default = "", type = str, help = "output directory for stratified gff files");
args = parser.parse_args();
attrs = args.attributes;

#parse options
if(args.output):
	dirname = args.output;
	try:
		os.mkdir(dirname)
	except:
		sys.stderr.write("directory %s already exists, gff files will be (over)written there\n" % os.path.abspath(dirname));		
else:
	dirname = args.output


stratified = {};
type2filehandler = {};


for interval in BedTool(args.path):
	type_ = tuple([interval.attrs[x] for x in attrs]);
	if(type_ in type2filehandler):
		stratified[type_]+=1;
	else:	
		filename = ".".join([os.path.splitext(os.path.basename(args.path))[0]] + list(type_) + ['gff'])
		type2filehandler[type_] = open(os.path.join(dirname, filename), 'w');		
		stratified[type_]=1;
	type2filehandler[type_].write(str(interval))
	

	
for handler in type2filehandler.values():
	handler.close();


sys.stderr.write("\t".join(args.attributes + ['total', 'fraction\n\n']));
total = float(sum(stratified.values()))
for name, count in stratified.items():
	sys.stderr.write("\t".join( list(name) + ["%d\t%1.4f\n" % (count, count/total)] ));

############################################################################################################################









