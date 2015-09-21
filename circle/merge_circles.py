#! /usr/bin/python
'''Merges circular reads into circles''' 
import sys;
import argparse
from collections import defaultdict, Counter

from pybedtools import BedTool

from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Merges circular reads into circles');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the circles, bed/gff file");
parser.add_argument('-od', '--dictionary', nargs = '?', required = True, type = str, help = "path to output \"interaction to read id\" file")
parser.add_argument('--support', nargs = '?', default = 2, type = int, help = "min read support for a circle to pass")
args = parser.parse_args()


def merge(reads, cnum):
	score = str(max([int(x.score) for x in reads]));
	attrs = [('n_uniq', len(reads)), ("ID", "c%d" % cnum)];	
	return construct_gff_interval(reads[0].chrom, reads[0].start, reads[0].stop, feature='circle', score=score, strand=reads[0].strand, source='.', frame='.', attrs=attrs)

circles = defaultdict(list);
for interval in BedTool(args.path):
	circles[(interval.chrom, interval.strand, interval.start, interval.stop)].append(interval);
	
cnum = 0;
passed = 0;
total = 0;
with open(args.dictionary, 'w') as fdict:
	for k, reads in circles.iteritems():
		total+=len(reads)
		if(len(reads)>=args.support):
			cnum+=1;
			passed+=len(reads);
			sys.stdout.write(str(merge(reads, cnum)));
			for read in reads:
				fdict.write("c%d\t%s\n" % (cnum, read.name));
				
				
				
sys.stderr.write("total circular reads: %d\npassed circular reads: %d\nfraction passed %1.5f\n\n" % (total, passed, float(passed)/total));
sys.stderr.write("%d circles generated from %d circular reads\n\n" % (cnum, passed));

sys.stderr.write("%s\n" % Counter([len(x) for x in circles.values()]).values());

