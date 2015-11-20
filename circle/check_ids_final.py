#! /usr/bin/python
'''Checks if reads which get rise to circles are actually in single bed circle file''' 
import argparse
import sys;
from collections import defaultdict, Counter

from Bio import SeqIO
from pybedtools import BedTool



parser = argparse.ArgumentParser(description='Checks if reads which get rise to circles are actually in single bed circle file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the chimeras, double bed/gff file");
parser.add_argument('--cid2rid', nargs = '?', required = True, type = str, help = "path to a dictionary file, which connects read ids with circular ids, chiflex tsv format");
parser.add_argument('--circle2read', nargs = '?', required = True, type = str, help = "path to a dictionary file, which connects circular ids with read ids, it is based on. find_circ.pt fasta format");
parser.add_argument('--circles', nargs = '?', required = True, type = str, help = "path to circular RNAs. bed/gff custom format");
parser.add_argument('--difference', nargs = '?', default = 0, type = int, help = "max allows distance between corresponding edges to count circles the");
args = parser.parse_args();

def compare(circle1, circle2, diff=8):
	
	start_diff = abs(circle1.start - circle2.start);
	stop_diff = abs(circle1.stop - circle2.stop);
	
	if(stop_diff<=diff and start_diff<=diff):
		return True
	else:
		return False


circles_fc = dict([(x.name, x) for x in BedTool(args.circles)])

rid2cid_fc = {};
rid2seq = {}
for seqrecord in SeqIO.parse(args.circle2read, 'fasta'):
	cid, readid = seqrecord.description.split(" ")
	if(cid in circles_fc):
		rid2cid_fc[readid] = cid;
		rid2seq[readid] = str(seqrecord.seq.upper())
	
	
circles_chiflex = dict([(x.name, x) for x in BedTool(args.path)])

rid2cid_chiflex = {};
with open(args.cid2rid) as f:
	for l in f:
		cid, readid = l.strip().split("\t")
		rid2cid_chiflex[readid] = cid


shared_reads = [];
shared_circles = set();
for rid, cid_fc in rid2cid_fc.iteritems():
	cid_chiflex = rid2cid_chiflex.get(rid, None);
	if(cid_chiflex):
		shared_reads.append(rid);
		intersection = compare(circles_fc[cid_fc], circles_chiflex[cid_chiflex], diff=args.difference)
		if(intersection):
			shared_circles.add(cid_fc);
		#else:
			#cfc = circles_fc[cid_fc]
			#ccf = circles_chiflex[cid_chiflex]
			
			#print "%s\t%d\t%d\t%s" % (cfc.chrom, cfc.start, cfc.end, cfc.strand)
			#print "%s\t%d\t%d\t%s\n" % (ccf.chrom, ccf.start, ccf.end, ccf.strand)
			#print rid2seq[rid]
			#print rid
			#print
			#print "*"*140
			
			#sys.stdout.write(str(circles_fc[cid_fc]))
			#sys.stdout.write("%s" % str(circles_chiflex[cid_chiflex]))
		
		
		
total_reads, found_reads = len(rid2cid_fc), float(len(shared_reads))
sys.stderr.write("total number of circular reads:\t%d\nnumber of circular reads detected:\t%d\nratio:\t%1.5f\n\n" % (total_reads, found_reads, float(found_reads)/total_reads));	
total_circles, found_circles = len(circles_fc), float(len(list(shared_circles)))
sys.stderr.write("total number of circles:\t%d\nnumber of circles detected:\t%d\nratio:\t%1.5f\n\n" % (total_circles, found_circles, float(found_circles)/total_circles));

sys.stderr.write("total number of chiflex circles:\t%d\nnumber of chiflex circular reads detected:\t%d\n\n" % (len(circles_chiflex), len(rid2cid_chiflex)));
