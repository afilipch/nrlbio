#! /usr/lib/python
'''Tool for destructive miRNA sites prediction'''

import argparse
import os
import sys
from collections import defaultdict

from Bio import SeqIO;
from Bio.Seq import reverse_complement;

from nrlbio.mirna import fasta2mirnas, MODES_ORDER



parser = argparse.ArgumentParser(description='Tool for destructive miRNA sites prediction');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the sequences to find binding sites in, fasta format");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the miRNAs, fasta format");
parser.add_argument('--slength', nargs = '?', default=7, type = int, help = "Length of the seed match, start position is always equals 2");
parser.add_argument('--alength', nargs = '?', default=20, type = int, help = "Length of the extension for miRNA target candidate");
parser.add_argument('--reverse', nargs = '?', default=False, const = True, type = bool, help = "If set, search will be performed also on reverse strand");
args = parser.parse_args();


match2mirid = defaultdict(list);
mirdict = {}

for mirid, mirna in fasta2mirnas(args.mir, seed_start=1, seed_stop=1+args.slength).items():
	mirdict[mirid] = mirna.seq;
	mirna.set_1mm_matches()
	match2mirid[mirna.match].append(mirid);
	for mm in mirna.matches_1mm:
		match2mirid[mm].append(mirid);
		
if(args.reverse):
	reverse2mirid = {};
	for match, mirids in match2mirid.items():
		reverse2mirid[reverse_complement(match)] = mirids

		
#start = 0;
end = args.slength

for seqrecord in SeqIO.parse(args.path, 'fasta'):
	l=[str(x) for x in seqrecord[:end]];
	mirids = match2mirid.get(''.join(l));
	if(mirids):
		print "%s\t%d\t%d\t%s\t0\t+" % (seqrecord.id, 0, end, ",".join(mirids))

	for c, n in enumerate(seqrecord[end:]):
		l.pop(0)
		l.append(n)
		mirids = match2mirid.get(''.join(l));
		if(mirids):
			print "%s\t%d\t%d\t%s\t0\t+" % (seqrecord.id, c+1, c+end+1, ",".join(mirids))

