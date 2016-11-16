#! /usr/bin/python
'''Checks if the chimeric read carry  information of nontemplate addition of nucleotides to miRNAs''' 
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool
import pysam;
from Bio import SeqIO



parser = argparse.ArgumentParser(description='Checks if the chimeric read carry  information of nontemplate addition of nucleotides to miRNAs');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "Path to the mapping hits for the first and second run");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the miRNA sequences, fasta file");
parser.add_argument('--targets', nargs = '?', required =True, type = str, help = "Path to the targets.gff (interactions in singlebed format) file from Chiflex. Statistics will be calculated only for the chimeras supporting these interactions");
parser.add_argument('--table', nargs = '?', required =True, type = str, help = "Path to the rid2iid.tsv file from Chiflex");
args = parser.parse_args();

### Get read ids 

iids = set();
for interval in BedTool(args.targets):
	iids.add(interval.attrs['ID']);
	
rids = set()
with open(args.table) as f:
	for l in f:
		rid, iid = l.strip().split('\t');
		if(iid in iids):
			rids.add(rid);


### Connect mappings from the first and the second rounds off mapping

mirid2len = {};
for seqrecord in SeqIO.parse(args.mir, 'fasta'):
	mirid2len[seqrecord.name] = len(seqrecord)


first = pysam.Samfile(args.path[0]);
second = pysam.Samfile(args.path[1]);

rid2pair = defaultdict(list);
rid2mirid = {};

for segment in first.fetch(until_eof=True):
	rid = segment.query_name 
	if((not rids) or rid in rids):
		rid2pair[rid].append(segment);
		rid2mirid[rid] = first.getrname(segment.tid);

for segment in second.fetch(until_eof=True):
	rid = segment.query_name 
	if((not rids) or rid in rids):
		rid2pair[rid].append(segment);
		
		
		
### Get statistics
		
cutdict = defaultdict(float)
addnucl = defaultdict(float)
addlen = defaultdict(float)		
		
		

for rid, (s1, s2) in rid2pair.iteritems():
	mirlen = mirid2len[rid2mirid[rid]]
	cut = mirlen - s1.reference_end
	cutdict[cut] += 1
	
	nucl = s2.query_sequence[:s2.query_alignment_start]
	addnucl[nucl] += 1
	addlen[len(nucl)] += 1
		
		
print addnucl;
print cutdict;
print addlen;
		
		
first.close();
second.close();