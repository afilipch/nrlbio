#! /usr/bin/python
'''Converts file from targetscan format to bed6 file'''
import sys;
import argparse;
import os

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts file from targetscan format to bed6 file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = os.path.abspath, help = "path to a directory with all target scan prdictions");
parser.add_argument('--specie', nargs = '?', required=True, type = str, help = "Name of the mirbase specie[hsa, mmu, cel and so on]");
args = parser.parse_args();


for path in [os.path.join(args.path, x) for x in os.listdir(args.path)]:
	for interval in BedTool(path):
		#print "-"*140
		#print interval
		processed = []
		interval.score = '0'
		gene, mirids = interval.name.split(":");
		mirids = mirids.split('/')
		processed.append('-'.join((args.specie, mirids[0].split('.')[0])))
		for mirid in mirids[1:]:
			processed.append('-'.join((args.specie, 'miR', mirid.split('.')[0])))
		for mirid in processed:
			interval.name = mirid
			print "\t".join(interval[:6])
		
	
		