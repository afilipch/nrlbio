# /usr/bin/python
'''converts genbank records into fasta files'''

import argparse
import sys;
import os;

from Bio import SeqIO

from nrlbio.ncbi import adjust_name



parser = argparse.ArgumentParser(description='converts genbank records into fasta files');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to genbank records");
parser.add_argument('--output', nargs = '?', required = True, type = str, help = "path to output folder");
args = parser.parse_args();

def change_names(path):
	for seq_record in SeqIO.parse(path, "genbank"):
		name = adjust_name(seq_record.name);
		seq_record.id = name
		seq_record.name = name
		seq_record.description = name
		yield seq_record



for path in args.path:
	bn = os.path.basename(path);
	sys.stderr.write("genbank file %s is being processed now\n" % bn)
	bn = ".".join(bn.split(".")[:-1] + ['fa'])
	fasta_name = os.path.join(os.path.abspath(args.output), bn)
	count = SeqIO.write(change_names(path), fasta_name, "fasta")
