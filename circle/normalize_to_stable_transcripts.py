#! /usr/bin/python
'''Normalizes an expression of circular splice junctions to the expression of linear splice junction of transcripts with a uniform expression(in TPMs) throughout different conditions(replicates)'''


import sys;
import argparse
from collections import defaultdict, namedtuple
from itertools import combinations

from pybedtools import BedTool
import numpy



parser = argparse.ArgumentParser(description='Normalizes an expression of circular splice junctions to the expression of linear splice junction of transcripts with uniform expression(in TPMs) throughout different conditions(replicates)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the cicrles, bed/gff format");
parser.add_argument('--lsj', nargs = '?', required = True, type = str, help = "Path to the linear splice junctions found in the same sample as given circles, bed/gff format");
parser.add_argument('--stable', nargs = '?', required = True, type = str, help = "Path to the linear splice junctions of transcripts with a uniform expression, bed/gff format");
#parser.add_argument('--maxvarcoeff', nargs = '?', default = 0.1, type = float, help = "Maximum variation coefficient allowed");
#parser.add_argument('--maxdev', nargs = '?', default = 0.2, type = float, help = "Maximum deviation from a mean allowed");
#parser.add_argument('--minexpr', nargs = '?', default = 50, type = float, help = "Min expression allowed");
args = parser.parse_args();

nice_factor = 1000000

lsj = BedTool(args.lsj);
stable = BedTool(args.stable);
norm_factor = 0

for interval in lsj.intersect(b = stable, f=0.99, F=0.99, u=True, s=True):
	#sys.stdout.write(str(interval));
	norm_factor += float(interval.attrs['n_uniq'])
	
sys.stderr.write("%1.1f\n" % norm_factor)	
	
norm_factor = norm_factor/nice_factor;
	
	
for interval in BedTool(args.path):
	interval.attrs['norm_expr'] = '%1.2f' % (float(interval.attrs['n_uniq'])/norm_factor)
	interval.attrs['Name'] = interval.name
	sys.stdout.write(str(interval))