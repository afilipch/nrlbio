# /usr/bin/python
'''adds features presented in gff/bed files to the genbank SeqRecord'''

import argparse
import sys;
from collections import defaultdict;

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio import SeqIO
from pybedtools import BedTool, Interval;



parser = argparse.ArgumentParser(description='converts genbank records into bed file usable for annotation');
parser.add_argument('-r', '--record', nargs = '?', type = str, help = "path to genbank record");
parser.add_argument('-f', '--features', nargs = '?', type = str, help = "path to the gff/bed file with features");
parser.add_argument('-o', '--output', nargs = '?', type = str, help = "path to the new genbank record with added features");
parser.add_argument('-ca', '--chr_adjustment', nargs = '?', default = '', type = str, help = "adds string to the name of chromosome in genbank file to find corresponding features in gff");
args = parser.parse_args();

sd = {'+': 1, '-': -1, '.': None};

def gff2feature(gff):
	location = FeatureLocation(gff.start, gff.end, strand = sd[gff.strand]);
	qualifiers = { "regulation_id": [gff.attrs['ID']], "note": [gff[2]] };
	return SeqFeature(location=location, type='regulation', qualifiers=qualifiers)


def add2seqrecord(gff_dict, genbank_file):
	for seq_record in SeqIO.parse(genbank_file, "genbank"):
		chrom = "".join([args.chr_adjustment, seq_record.name]);
		for gff in gff_dict[chrom]:
			seq_record.features.append(gff2feature(gff));
		yield seq_record
		
		
gff_dict = defaultdict(list)
for gff in BedTool(args.features):
	gff_dict[gff.chrom].append(gff);

SeqIO.write(add2seqrecord(gff_dict, args.record), args.output, "genbank")

