#! /usr/bin/python	
'''intellectually merges interaction files'''
import interaction_lib;
import mir_lib;
import sys;
import os;
import argparse;
import random;
from collections import *;


parser = argparse.ArgumentParser(description='intellectually merges interaction files');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-r', '--read2int', nargs = '+', default = [], type = str, help = "path to reads to interaction mappings (output/read2int.tsv)");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
parser.add_argument('-n', '--name', nargs = '?', default = "i", type = str, help = "name of dataset")
args = parser.parse_args();

interactions = [];
for path in args.path:
	interactions += interaction_lib.get(path, True);

int2read = defaultdict(set);
for p in args.read2int:
	f = open(p);
	for l in f:
		a = l.strip().split("\t");
		int2read[a[1]].add(a[0])
	f.close();	

	
exec("from sequence_data.systems import %s as gsys" % args.system)	
merged, int2new = interaction_lib.merge(interactions, gsys, int2read, args.name);

for el in merged:
	print "\t".join([str(x) for x in list(el)]);
	
if(args.read2int):	
	o = open(args.name + "_read2int.tsv", 'w');
	for k, v in int2new.iteritems():
		for r in v:
			o.write("%s\t%s\n" % (r,k));
	o.close();
