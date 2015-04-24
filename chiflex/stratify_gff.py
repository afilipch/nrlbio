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
parser.add_argument('-o', '--output', nargs = '?', default = "", type = str, help = "output directory for strtified gff files");
args = parser.parse_args();

#parse options
if(args.output):
	dirname = args.output;
	try:
		os.mkdir(dirname)
	except:
		sys.stderr.write("directory %s already exists, gff files will be (over)written there\n" % os.path.abspath(dirname));		
else:
	dirname = args.output


stratified = defaultdict(list);
total = 0;

for interval in BedTool(args.path):
	stratified[tuple([interval.attrs[x] for x in args.attributes])].append(interval);
	total += 1
	
str_int = dict([(x, len(stratified[x])) for x in stratified])	
	
for name, sf in stratified.iteritems():
	filename = ".".join([os.path.splitext(os.path.basename(args.path))[0]] + list(name) + ['gff'])
	with open(os.path.join(dirname, filename), 'w') as f:
		for interval in sf:
			f.write(str(interval))


sys.stderr.write("\t".join(args.attributes + ['total', 'fraction\n\n']));
for name, count in str_int.iteritems():
	sys.stderr.write("\t".join( list(name) + ["%d\t%1.4f\n" % (count, count/float(total))] ));

############################################################################################################################









