#! /usr/lib/python
'''orders bed file to doublebed format, that is connected pairs of intervals are consecutive''' 
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool

#a = [[1, 7], [4, 6], [1, 3], [4, 0], [0, 0]]
#for i in groupby(sorted(a, key=itemgetter(0)), key=itemgetter(0)):
	#for x in i[1]:
		#print i[0], x;
#sys.exit()

parser = argparse.ArgumentParser(description='orders bed file to doublebed format, that is: connected pairs of intervals are consecutive in the ordered file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file to be annotated");
parser.add_argument('--convert',nargs = '?', default = False, const = True, type = bool, help = "converts input file itself");
args = parser.parse_args();


pairs = defaultdict(list);
for interval in BedTool(args.path):
	pairs[interval.name.split("|")[0]].append(interval);

if(args.convert):
	with open(args.path, 'w') as f:
		for name, pair in sorted(pairs.items(), key = lambda x: int(x[0].split("_")[-1])):
			for interval in pair:
				f.write(str(interval));
else:		
	for name, pair in sorted(pairs.items(), key = lambda x: int(x[0].split("_")[-1])):
		for interval in pair:
			sys.stdout.write(str(interval));
