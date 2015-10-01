#! /usr/lib/python
'''Merges mapping hits into clusters(peaks)'''
import argparse
import sys;
import os;

from pybedtools import BedTool;

from nrlbio.mapcluster import Cluster;


parser = argparse.ArgumentParser(description='Merges mapping hits into clusters(peaks)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the bedgraph file");
parser.add_argument('-n', '--support', nargs = '?', default = 2, type = int, help = "min coverage support for the cluster");
parser.add_argument('-d', '--drop', nargs = '?', default = 2.5, type = float, help = "min drop of coverage to stop cluster extension");
parser.add_argument('-s', '--strand', nargs = '?', default = '.', type = str, help = "strand of the clusters");
args = parser.parse_args();
support = args.support;
drop = args.drop;
strand = args.strand;

startnew=True;
count = 1;
for interval in BedTool(args.path):
	score = int(interval[3]);
	if(startnew):
		if(score>=support):
			start = interval.start;
			chrom = interval.chrom;
			prevcov = score;
			end = interval.end;
			startnew = False;
			count += 1;
		else:
			pass;
	else:
		if(score<support or interval.start!=end or chrom!=interval.chrom or prevcov>score*drop):
			sys.stdout.write(str(Cluster("mc%d" % count, chrom, start, end, strand).gff()));
			if(score>=support):
				start = interval.start;
				chrom = interval.chrom;
				prevcov = score;
				end = interval.end;
				count += 1;
			else:
				startnew=True;
		else:
			prevcov = score;
			end = interval.end;

			
		







#parser = argparse.ArgumentParser(description='Merges mapping hits into clusters(peaks)');
#parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the sam/bam  file");
#parser.add_argument('-n', '--support', nargs = '?', default = 2, type = int, help = "min coverage support for the cluster");
#parser.add_argument('-d', '--drop', nargs = '?', default = 2.5, type = float, help = "min drop of coverage to stop cluster extension");
##parser.add_argument('-o', '--output', nargs = '?', default = "", type = str, help = "output directory for stratified gff files");
#args = parser.parse_args();
#support = args.support;
#drop = args.drop;


#samfile = pysam.AlignmentFile(args.path)
#startnew = True;
#segments = [];
#count=0;
#for c, pileupcolumn in enumerate(samfile.pileup()):
	#if(startnew):
		#if(pileupcolumn.n>=support):
			#start=pileupcolumn.pos;
			#chrom = samfile.getrname(pileupcolumn.reference_id)
			#prevcov=pileupcolumn.n;
			#count+=1;
			#startnew=False;
		#else:
			#pass;
	#else:
		#if(pileupcolumn.n<support or prevcov>pileupcolumn.n*args.drop):
			#end = pileupcolumn.pos;
			#startnew = True;
			#cluster = Cluster("mc%d" % count, chrom, start, end, segments)
			#sys.stdout.write(str(cluster.gff()));
		#else:
			#prevcov=pileupcolumn.n;
			
			
	#if(c>1000000):
		#break;
		
		
		
		
	

