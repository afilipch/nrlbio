#! /usr/lib/python
'''Generates markov models from provided sequences in fasta format'''

import argparse
import os
import sys

from Bio import SeqIO

from nrlbio.HMM import MultiMarkov, CONF_MULTIMARKOV

parser = argparse.ArgumentParser(description='Generates markov models from provided sequences in fasta format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file");
parser.add_argument('--order', nargs = '?', default = 1, type = int, help = "Markov model order. Correspondes to shuffling of \'order\'-nucleotides");
parser.add_argument('--output', nargs = '?', required=True, type = str, help = "Path to the folder, where models will be stored as yml files. The directory will be created, if doesn't exist");
args = parser.parse_args();




directory=os.path.abspath(args.output)
if not os.path.exists(directory):
    os.makedirs(directory)
    sys.stderr.write("directory %s is created, files will be written there\n\n" % directory)
else:
	sys.stderr.write("%s is already exists, files will be written there\n\n" % directory);

for seqrecord in SeqIO.parse(args.path, "fasta"):
	mm = MultiMarkov.from_seqrecord(seqrecord, args.order, **CONF_MULTIMARKOV)
	mm.serialize(os.path.join(directory, "%s.yml" % seqrecord.name))
sys.exit();