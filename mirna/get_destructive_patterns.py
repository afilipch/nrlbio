#! /usr/bin/python
'''Selects for targets with an extensive mirna sequence complementarity''' 
import argparse
import sys;
import math;

from nrlbio.generators import generator_doublebed



parser = argparse.ArgumentParser(description='Outputs targets of miRNA');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the mirna:target interactions to select from, double gff/bed. Pattern and energy have to be assigned to each interaction");
args = parser.parse_args();

def pattern2score(pattern):
	p = [int(x) for x in pattern.split(",")]
	return sum(p[1:9]) - sum(p[9:11]) + sum(p[11:19])

def interval2score(interval):
	return pattern2score(interval.attrs['pattern']) - float(interval.attrs['energy']) + math.log(int(interval.attrs['n_uniq']));

	
for i1, i2 in generator_doublebed(args.path):
	i2.name = i1.chrom
	i2.attrs['dscore'] = str(interval2score(i2))
	sys.stdout.write(str(i2))