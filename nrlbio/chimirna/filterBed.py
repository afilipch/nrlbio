#! /usr/bin/python	
'''Script filter bed file according to a present of seed complementarity to a miRNA in provided file'''
#import sam_lib;
import mir_lib
import gconst;
import chipart_lib;
import cluster_lib
import parsing;
import sys;
import os;
import argparse;
import logging;
from collections import *;



parser = argparse.ArgumentParser(description='Script filter bed file according to a present of seed complementarity to a miRNA in provided file');
# input files
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to the bed-like or interactions like file");
parser.add_argument('--mir', nargs = '?', type = str, help = "filter reads which have seed match to miRNA provided in the list");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
parser.add_argument('-i', '--interactions', nargs = '?', default = False, const = True, type = str, help = "genome system")
args = parser.parse_args();

Bed = namedtuple("Bed","chromosome, start, end, bid, score, strand")
exec("from sequence_data.systems import %s as gsys" % args.system);

mirdict = parsing.fasta2dict(gconst.system2mir[args.system]);
if(args.mir):
	mirselect = [];
	f = open(args.mir)
	for l in f:
		mirselect.append(l.strip());
	f.close()	
	mirdict =  dict([(k, mirdict[k]) for k in filter(lambda x: x in mirselect, mirdict.keys())])	

seeds = mir_lib.seed2mirs(mirdict)
mseeds = mir_lib.mir2seed(mirdict)


b = open(args.path[0]);
for l in b:
	a = l.strip().split("\t");
	bed = Bed(a[0], int(a[1]), int(a[2]), a[3], int(a[4]), a[5])
	seq = gsys.genome.get_oriented(bed.chromosome, bed.start, bed.end, bed.strand).upper()
	if(args.interactions):
		intersection = set(a[6].split(",")) & set(mirselect)
		if(intersection):
			mirid = list(intersection)[0];
			if(True or mseeds[mirid] in seq):
				print l.strip()
			
	else:	
		for k,v in seeds.iteritems():
			if k in seq:
				print l.strip()
				break;
b.close()			



