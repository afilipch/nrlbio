#! /usr/bin/python
'''Provides exonic structure of given circles on basis of transcript annotation(gff3 file) in bed12 format'''


import sys;
import argparse
from collections import defaultdict
from itertools import combinations

from pybedtools import BedTool

from nrlbio.itertools_extension import powerset



parser = argparse.ArgumentParser(description='Provides exonic structure of given circles on basis of transcript annotation(gff3 file) in bed12 format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the exons of circles, bed format");
args = parser.parse_args();



def check_consistency(exons, gap):
	if(len(exons)<2):
		return True
	else:
		for e1, e2 in zip(exons, exons[1:]):
			if(e1.end + gap > e2.start):
				return False
		else:
			return True
		
		
def get_route(exons, gap):
	routes = defaultdict(list)
	for p1, p2 in combinations(range(len(exons)), 2):
		if((exons[p1].end + gap) < exons[p2].start):
			routes[p2].append([p1])
			for el in routes[p1]:
				routes[p2].append(el+[p1])
	
	ans = [[x] for x in exons]
	for k, l in routes.iteritems():
		for v in l:
			ans.append(v+[k])
			
	return ans;
		
	

		 
		 
		 


def exons2bed12(exons, gap = 20):
	start = exons[0].start
	end = exons[-1].end
	start_exons = [x for x in exons if x.start==start]
	end_exons = [x for x in exons if x.end==end]
	inner_exons = [x for x in exons if x.start!=start and x.end!=end]
	
	if(len(inner_exons)<30):
		get_route(inner_exons, gap)
		print 1;	

	
	
circles = defaultdict(list)
for interval in BedTool(args.path):
	if(interval.score == '1'):
		circles[interval.name].append(interval);
		
		

for name, exons in circles.iteritems():
	exons2bed12(exons)
