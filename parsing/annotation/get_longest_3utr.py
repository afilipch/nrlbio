#! /usr/lib/python
'''extracts longest UTR for each coding gene in genbank records''' 
import argparse
import sys;
from collections import defaultdict;

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio import SeqIO

from nrlbio.ncbi import flanking_exons, adjust_name






parser = argparse.ArgumentParser(description='extracts longest UTR for each coding gene in genbank records');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to genbank records");
parser.add_argument('-t', '--type', nargs = '?', choices = ['utr3', 'utr5'], required = True, type = str, help = "type of utr");
parser.add_argument('-f', '--format', nargs = '?', choices = ['fasta', 'bed'], required = True, type = str, help = "format of output");
args = parser.parse_args();

sd = {1: '+', -1: '-', None: '.'};



def split_mrna(mrna, cds):
	strand = mrna.location.strand
	
	if(type(mrna.location) == CompoundLocation):
		exons = mrna.location.parts;
	else:
		exons = [mrna.location];
	
	left_exons, right_exons = flanking_exons(exons, cds.location.start, cds.location.end, strand)

	# we always return exons in an order 5'UTR, 3'UTR;		
	if(strand == 1):
		return left_exons, right_exons
	elif(strand == -1):
		return right_exons, left_exons
		
		
		

def compound_to_bed(exons, chrom, name):
	for e in exons:
		print "%s\t%d\t%d\t%s\t0\t%s" % (chrom, e.start, e.end, name, sd[e.strand]);
	return 1	
	
		
		
def fetch_utr(seq_record):
	gene2isoforms = defaultdict(list);
		
	for feature in seq_record.features:
		
		if(feature.type == "mRNA"):
			mrna = feature;
		elif(feature.type == "CDS"):
			utr5, utr3 = split_mrna(mrna, feature)
			if(args.type == 'utr3' and utr3):
				gene2isoforms[mrna.qualifiers['gene'][0]].append(utr3);
			elif(args.type == 'utr5' and utr5):
				gene2isoforms[mrna.qualifiers['gene'][0]].append(utr5);
			else:
				pass;
		else:
			pass;
			
	return gene2isoforms;
	
	
def get_longest(seq_record, gene2isoforms):
	l = []
	c = 0;

	chrom = adjust_name(seq_record.name);	
	
	for gene, isoforms in gene2isoforms.iteritems():
		longest = max(isoforms, key = lambda i: sum([len(x) for x in i]))
		
		if(args.format == 'bed'):
			compound_to_bed(longest, chrom, gene)
			
		elif(args.format == 'fasta'):
			if(len(longest) > 1):
				location = CompoundLocation(longest, operator = "join")
			else:
				location = longest[0];
				
			feature = SeqFeature(location=location, type='utr', strand = longest[0].strand)
			
			#print longest[0].strand
			
			f = feature.extract(seq_record)
			f.name = gene
			f.id = gene
			f.description = gene
			l.append(f);
	return l;	
		

def run(path):	
	for p in path:
		for seq_record in SeqIO.parse(p, "genbank"):
			sys.stderr.write("sequence record %s is being processed now\n" % seq_record.name)
			gene2isoforms = fetch_utr(seq_record);	
			for f in get_longest(seq_record, gene2isoforms):
				yield f;
			
			
if(args.format == 'bed'):
	run(args.path);
elif(args.format == 'fasta'):
	SeqIO.write(run(args.path), sys.stdout, 'fasta')
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			