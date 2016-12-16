#! /usr/lib/python
'''assigns sequence to each gff entry'''

import argparse
import sys

from Bio import SeqIO;
from pybedtools import BedTool

from nrlbio.pybedtools_extension import interval2seq

parser = argparse.ArgumentParser(description='assigns sequence to each gff entry');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to gff file to be annotated");
parser.add_argument('-f', '--fasta', nargs = '+', required = True, type = str, help = "Path to reference(genome, mirBase, transcriptome) fasta files");
args = parser.parse_args();


def seqrecord_iterator(paths):
	for path in paths:
		for seqrecord in SeqIO.parse(path, "fasta"):
			yield seqrecord;


reference = SeqIO.to_dict(seqrecord_iterator(args.fasta))

#sys.stderr.write("\n".join([str(x) for x in reference.keys() if x[:3]=='chr']) + "\n")

for interval in BedTool(args.path):
	seq = interval2seq(interval, reference);
	interval.attrs['seq'] = seq;
	sys.stdout.write(str(interval))
	