#! /usr/bin/python	
'''Script filter miRNA fasta file to top expressed entries'''
#import sam_lib;
import rnahybrid;
import argparse;
from Bio import SeqIO;



parser = argparse.ArgumentParser(description='Script filter miRNA fasta file to top expressed entries');
# input files
parser.add_argument('-e', '--expressed', nargs = '?', type = str, help = "custom file of mirna expression");
parser.add_argument('--mir', nargs = '?', type = str, help = "miRNA for organism of interest as fasta file");
parser.add_argument('-n', '--top', nargs = '?', default = 100, type = int, help = "show n top expressed miRNA");
args = parser.parse_args();


mirdict = SeqIO.to_dict(SeqIO.parse(args.mir, "fasta"))
expressed = [];

f = open(args.expressed);
for l in f:
	a = l.strip().split();
	expressed.append((a[1], int(a[0])));
f.close()	
expressed.sort(key = lambda x: x[1], reverse = True)
myfilter = [x[0] for x in expressed[:args.top]];

for mirid in myfilter:
	print ">%s" % mirid
	print mirdict[mirid].seq


