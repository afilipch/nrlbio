#! /usr/bin/python
'''Checks if mappings to small rna reference carry information of nontemplate addition of nucleotides to smallRNAs''' 
import argparse
import sys;
import os;
from collections import defaultdict, Counter

import pysam;
from Bio import SeqIO



parser = argparse.ArgumentParser(description='Checks if mappings to small rna reference carry information of nontemplate addition of nucleotides to smallRNAs');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the mapping hits to small rna reference");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the downstream-extended miRNA sequences, fasta file");
parser.add_argument('--mirids', nargs = '+', default=None, type = str, help = "If set, only the provided miRNAs will be analyzed");
parser.add_argument('--output', nargs = '?', required=True, type = str, help = "Path to the output folder");
args = parser.parse_args();

###miRNA processing
mirid2seq = defaultdict(set)
for seqrecord in SeqIO.parse(args.mir, 'fasta'):
	mirid2seq[seqrecord.name].add(str(seqrecord.seq.upper()))
	
if(args.mirids):
	selected_mirids = args.mirids;
else:
	selected_mirids = mirid2seq.keys()	


###Samfile processing
mir_pos_to_nucls = defaultdict(lambda: defaultdict(lambda: defaultdict(float)));

for path in args.path:
	samfile = pysam.Samfile(path);
	for segment in samfile.fetch(until_eof=True):
		mirid = samfile.getrname(segment.tid);
		start = segment.reference_start;
		
		right = segment.query_sequence[segment.query_alignment_start:]
		for pos, nucl in enumerate(right):
			mir_pos_to_nucls[mirid][pos+start][nucl] += 1;
		
		left = segment.query_sequence[max(0, segment.query_alignment_start - start) :segment.query_alignment_start]
		for pos, nucl in enumerate(left[::-1]):
			mir_pos_to_nucls[mirid][start-pos][nucl] += 1;
	samfile.close();

###outputting
norder = 'A', 'C', 'T', 'G';


for mirid in selected_mirids:
	with open(os.path.join(args.output, mirid + '.tsv'), 'w') as f:
		f.write("position\t%s\tmiRNA_precursor\n" % "\t".join(norder))
		mirdict = mir_pos_to_nucls[mirid];
		mirseqs = list(mirid2seq[mirid]);
		porder = range(max([len(x) for x in mirseqs]))

		for pos in porder: 
			d = mirdict[pos];
			if(not d):
				d = {'A':0, 'C':0, 'T':0, 'G':0}
				
			precnucls = ",".join([x[pos] for x in mirseqs])
			f.write("%d\t%s\t%s\n" % (pos+1, "\t".join(["%1.1f" % d[x] for x in norder]), ",".join([x[pos] for x in mirseqs]) ) )
		
		
