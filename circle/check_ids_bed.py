#! /usr/bin/python
'''Checks if reads which get rise to circles are actually in double bed/gff file''' 
import argparse
import sys;
from collections import defaultdict, Counter

from Bio import SeqIO
from pybedtools import BedTool



parser = argparse.ArgumentParser(description='Checks if reads which get rise to circles are actually in double bed/gff file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the chimeras, double bed/gff file");
parser.add_argument('--circle2read', nargs = '?', required = True, type = str, help = "path to a dictionary file, which connects circular ids with read ids, it is based on. fasta format");
parser.add_argument('--circles', nargs = '?', required = True, type = str, help = "path to circular RNAs. bed/gff custom format");
args = parser.parse_args();

def compare(circle, intervals, diff=8):
	if(not intervals):
		return False
	
	start_diff = min([abs(circle.start-x.start) for x in intervals]);
	stop_diff = min([abs(circle.stop-x.stop) for x in intervals]);
	
	if(stop_diff<=diff and start_diff<=diff):
		return True
	else:
		return False


circles = dict([(x.name, x) for x in BedTool(args.circles)])


circid2read = defaultdict(list);
for seqrecord in SeqIO.parse(args.circle2read, 'fasta'):
	cid, readid = seqrecord.description.split(" ")
	circid2read[cid].append(readid);
	
	
	

readid2intervals = defaultdict(list)
for l in circid2read.values():
	for readid in l:
		readid2intervals[readid] = [];



for interval in BedTool(args.path):
	name = interval.name.split("|")[0]
	if(name in readid2intervals):
		readid2intervals[name].append(interval);




correct_circle = 0;
correct_read = 0;
total_filtered_read = 0;
for cid, circle in circles.iteritems():
	correct = 0;
	for readid in circid2read[cid]:
		correct += int(compare(circle, readid2intervals[readid]))
		total_filtered_read += 1;
	correct_circle += min(1, correct);
	correct_read += correct	

		
		
total_read = sum([len(x) for x in circid2read.values()]);
found_read = len(filter(len, readid2intervals.values()));
sys.stderr.write("total number of circular reads:\t%d\nnumber of circular reads detected:\t%d\nratio:\t%1.5f\n\n" % (total_read, found_read, float(found_read)/total_read));
sys.stderr.write("total number of circular reads:\t%d\nnumber of circular reads detected correctly:\t%d\nratio:\t%1.5f\n\n" % (total_filtered_read, correct_read, float(correct_read)/total_filtered_read));
sys.stderr.write("total number of circles:\t%d\nnumber of circles detected correctly:\t%d\nratio:\t%1.5f\n\n" % (len(circles), correct_circle, float(correct_circle)/len(circles)));
sys.stderr.write("num of hits\tnum of reads\n")
for kv in Counter([len(x) for x in readid2intervals.values()]).iteritems():
	sys.stderr.write("%d\t%d\n" % kv);
