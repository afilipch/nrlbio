#! /usr/bin/python
'''Format fasta file produced by circbase in a way compatible with chiflex''' 
import sys;
import argparse
from collections import defaultdict, Counter;


from Bio import SeqIO



parser = argparse.ArgumentParser(description='Format fasta file produced by circbase in a way compatible with chiflex');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the file with circles\' sequences produced by circbase(only exons included), fasta format");
parser.add_argument('--bed', nargs = '?', required = True, type = str, help = "path to the corresponding bed file")
parser.add_argument('--length', nargs = '?', default = 40, type = int, help = "min length of circle")
args = parser.parse_args();

def description2interval(description):
	name, nreads, rlength, rcoord, rstrand, gid, gsymbol, annotation = description.split(" ");
	chrom = rcoord.split(":")[0]
	start, end = rcoord.split(":")[1].split("-")
	score = nreads.split("=")[1]
	strand = rstrand.split("=")[1]
	length = rlength.split("=")[1]
	
	return "\t".join((chrom, start, end, name, score, strand, length, gid, gsymbol, annotation))


removed=0;
passed=0;

with open(args.bed, 'w') as f:
	for seqrecord in SeqIO.parse(args.path, "fasta"):
		if(len(seqrecord.seq)>=args.length and ('N' not in seqrecord.seq)):
			f.write("%s\n" % description2interval(seqrecord.description))
			print ">%s" % seqrecord.name
			print seqrecord.seq.upper();
			passed+=1
		else:
			removed+=1;
				

sys.stderr.write("total circles: %d\npassed circles: %d\n removed circles: %d\n" % (removed+passed, passed, removed))