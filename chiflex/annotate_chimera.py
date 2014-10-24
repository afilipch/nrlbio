#! /usr/lib/python
'''assignes the type(linera/circular splice junctions, intra/inter-molecular interaction) to each chimera/interaction.''' 
import argparse
import os;

import pysam;

from pybedtools import BedTool, Interval;



parser = argparse.ArgumentParser(description='assignes the type(linera/circular splice junctions, intra/inter-molecular interaction) to each chimera/interaction');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-e', '--exons', nargs = '?', default = "sam", type = str, help = "path to the bed file of transcripts\' exons");
args = parser.parse_args();

def list2interval(l):
	return Interval(l[0], int(l[1]), int(l[2]), l[3], l[4], l[5])

def chimera2intervals(path, n):
	with open(path) as f:
		for l in f:
			a=l.strip().split("\t");
			yield list2interval(a[6*n:6*(n+1)])


class annotated_chimera(object):
	def __init__(self, left_interval, right_interval, left_intersection, right_intersection):
		self.left_interval = left_interval
		self.right_interval = right_interval
		self.left_intersection = left_intersection
		self.right_intersection = right_intersection
		


			

def generate_linked_intervals(intersection1, intersection2):
	curname = '';
	curleft = None;
	left_intersection = [];
	right_intersection = [];
	
	for i in intersection1:
		if(i.name == curname):
			left_intersection.append(list2interval(i[6:12]));
		elif(curleft):
			
first_intervals = BedTool(chimera2intervals(args.path, 0))
second_intervals = BedTool(chimera2intervals(args.path, 1))
exons = BedTool(args.exons)

print second_intervals.intersect(exons, wao=True, s=True)

#intersect(b, wao=True)
#sg_0_01_124031_x5
#