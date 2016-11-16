#! /usr/bin/python
'''Checks if mappings to small rna reference carry information of nontemplate addition of nucleotides to smallRNAs''' 
import argparse
import sys;
from collections import defaultdict, Counter

import pysam;
from Bio import SeqIO



parser = argparse.ArgumentParser(description='Checks if mappings to small rna reference carry information of nontemplate addition of nucleotides to smallRNAs');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mapping hits to small rna reference");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the miRNA sequences, fasta file");
args = parser.parse_args();

#get miRNA precursor extension nucleotides 

### Connect mappings from the first and the second rounds off mapping

mirid2extension = defaultdict(set)
for seqrecord in SeqIO.parse(args.mir, 'fasta'):
	mirid2extension[seqrecord.name].add(str(seqrecord.seq.upper()[-5:]))


mirid_to_nucls = defaultdict(lambda: defaultdict(float));
mirid_to_len = defaultdict(lambda: defaultdict(float));

samfile = pysam.Samfile(args.path);
for segment in samfile.fetch(until_eof=True):
	mirid = samfile.getrname(segment.tid);
	count = float(segment.query_name.split('_')[-1][1:])
	right = segment.query_sequence[segment.query_alignment_end:]
	
	mirid_to_len[mirid][len(right)] += count
	
	if(right):
		count = count/len(right);
		for nucl in right:
			mirid_to_nucls[mirid][nucl] += count;
	else:
		mirid_to_nucls[mirid]['N'] += count;
			

mirid_to_mean_len = {};
for mirid, d in mirid_to_len.items():
	total = float(sum(d.values()))
	length = sum([x[0]*x[1] for x in d.items()])
	mirid_to_mean_len[mirid] = length/total

nucl_list = [];	
for mirid, ndict in mirid_to_nucls.items():
	total = sum(ndict.values())
	if(total > 100):
		nucl_list.append(tuple([mirid, total, mirid_to_mean_len[mirid]] + [ndict[x]/total for x in 'NACTG'] + [", ".join(mirid2extension[mirid])] ))
		
		
nucl_list.sort(key= lambda x: x[1], reverse=True);
print "\t".join(['miRNA', 'total reads', 'mean length of 3\' extension', 'fraction of no-extended, %', 'fraction of added "A", %', 'fraction of added "C", %', 'fraction of added "U", %', 'fraction of added "G", %', '3\' extension in precursors'])
for l in nucl_list:
	print "%s\t%d\t%1.1f\t%s\t%s" % (l[0], l[1], l[2], "\t".join([ "%1.1f" % (x*100) for x in l[3:8]]), l[8])

		
#print mirid_to_nucls['hsa-miR-671-5p']