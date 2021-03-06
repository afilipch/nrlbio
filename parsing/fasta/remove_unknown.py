# /usr/bin/python
'''removes entries of fasta file with uncalled (N) bases'''

import argparse
import sys;
from collections import defaultdict, Counter;

from Bio import SeqIO;

from nrlbio.itertools_extension import median;


parser = argparse.ArgumentParser(description='removes entries of fasta file with uncalled (N) bases');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file");
args = parser.parse_args();


def run(path):
	for seq_record in SeqIO.parse(path, 'fasta'):
		if('N' not in seq_record.seq.upper()):
			print ">%s" % seq_record.id
			print seq_record.seq.upper()
		else:
			pass
		
run(args.path)		
#SeqIO.write(run(args.path), sys.stdout, "genbank");
	


