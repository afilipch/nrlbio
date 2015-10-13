#! /usr/bin/python
'''Creates bed file of randomly chosen intervals from intervals of the parent bed file. The randomly chosen intervals length distribution has to be of the same size as in child bed file. Script is useful to check if the child intervals were chosen from parent ones randomly or not.'''
import argparse
import sys;
import os;
import random

from pybedtools import BedTool;
from nrlbio.random_extension import weighted_choice


parser = argparse.ArgumentParser(description='Creates bed file of randomly chosen intervals from intervals of the parent bed file. The randomly chosen intervals length distribution has to be of the same size as in child bed file. Script is useful to check if the child intervals were chosen from parent ones randomly or not.');
parser.add_argument('-p', '--parent', nargs = '?', required = True, type = str, help = "path to the parent intervals, bed/gff file");
parser.add_argument('-c', '--child', nargs = '?', required = True, type = str, help = "path to the child intervals, bed/gff file");
args = parser.parse_args();

parent_list = list(BedTool(args.parent))
ntrials = 1000;


#parent_dict = dict([(x.name, x) for x in parent_list]);

#def get_interval(interval, parent_dict):
	#parent_length = [(x[0], len(x[1])) for x in parent_dict.items() if len(x[1])>=len(interval)]
	#if(parent_length):
		#name = weighted_choice(parent_length)
		#pi = parent_dict[name];
		#start = pi.start + random.randint(0, len(pi) - len(interval))
		#return pi.chrom, str(start), str(start+len(interval)), interval.name, '0', pi.strand
	#else:
		#sys.stderr.write("Warning length of the child interval is greater than any of parent ones. The interval will be skipped\n")
		#return None;


def get_fast(interval, parent_list):
	for _ in range(ntrials):
		pi = random.choice(parent_list);
		if(len(pi)>len(interval)):
			start = pi.start + random.randint(0, len(pi) - len(interval))
			return pi.chrom, str(start), str(start+len(interval)), interval.name, '0', pi.strand		
	else:
		sys.stderr.write("Warning length of the child interval is greater than any of randomly chosen sample(%d) of parent ones. The interval will be skipped\n" % ntrials)
		return None



for ci in  BedTool(args.child):
	ri = get_fast(ci, parent_list)
	if(ri):
		print "\t".join(ri);
