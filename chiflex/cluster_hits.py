#! /usr/bin/python
'''Merges mapping hits into clusters(peaks)'''
import argparse
import sys;
import os;

from pybedtools import BedTool;
import pysam;

from nrlbio.mapcluster import Cluster;


parser = argparse.ArgumentParser(description='Merges mapping hits into clusters(mapping peaks)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the bedgraph file");
parser.add_argument('-b', '--sbam', nargs = '?', required = True, type = str, help = "path to the sorted sam/bam file");
parser.add_argument('-n', '--support', nargs = '?', default = 2, type = int, help = "min coverage support for the cluster");
parser.add_argument('-d', '--drop', nargs = '?', default = 2.5, type = float, help = "min drop of coverage to stop cluster extension");
parser.add_argument('-s', '--strand', nargs = '?', default = '.', choices = ['+', '-', '.'], type = str, help = "strand of the clusters");
args = parser.parse_args();

SUPPORT = args.support;
DROP = args.drop;
STRAND = args.strand;

def strand_callback(segment):
	if(STRAND == '+'):
		return not segment.is_reverse
	elif(STRAND == '-'):
		return segment.is_reverse
	else:
		return True;
	
	
startnew=True;
count = 1;
peak = 0;
raw_clusters = [];
for interval in BedTool(args.path):
	score = int(interval[3]);
	if(startnew):
		if(score>=SUPPORT):
			start = interval.start;
			chrom = interval.chrom;
			prevcov = score;
			end = interval.end;
			startnew = False;
			peak = score;
			count += 1;
		else:
			pass;
	else:
		if(score<SUPPORT or interval.start!=end or chrom!=interval.chrom or prevcov>score*DROP):
			raw_clusters.append(Cluster("mc%d" % count, chrom, start, end, peak, STRAND));
			if(score>=SUPPORT):
				start = interval.start;
				chrom = interval.chrom;
				prevcov = score;
				end = interval.end;
				peak = score;
				count += 1;
			else:
				startnew=True;
		else:
			prevcov = score;
			peak = max(peak, score);
			end = interval.end;
else:
	if(not startnew):
		raw_clusters.append(Cluster("mc%d" % count, chrom, start, end, peak, STRAND));


samfile = pysam.AlignmentFile(args.sbam);

for cluster in raw_clusters:
	for segment in samfile.fetch(reference=cluster.chrom, start=cluster.start, end=cluster.end):
		if(strand_callback(segment) and  (min(cluster.end, segment.reference_end) - max(cluster.start, segment.reference_start) > segment.reference_length*0.5)):
			cluster.segments.append(segment);
	cluster.get_support();		
	if(cluster.support>=SUPPORT):
		sys.stdout.write(str(cluster.gff()));
	

		
		
		
	

