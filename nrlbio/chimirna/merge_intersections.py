import random;
import sys;
import re;
import os;
import argparse;
import multiprocessing;
from collections import defaultdict;




parser = argparse.ArgumentParser(description='merge intersection with known interactions for one specie');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to intersections with known interactions file");
args = parser.parse_args();

d = defaultdict(set);
c = 0;
for path in args.path:
	f = open(path);
	for l in f:
		c += 1;
		a = l.strip().split("\t");
		d[tuple(a[:2])].add(a[2])
	f.close();
	
for k in sorted(d.keys()):
	v = d[k]
	print "\t".join((k[0], k[1], ",".join(v)))
#print c