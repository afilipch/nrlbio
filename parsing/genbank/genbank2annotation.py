#! /usr/lib/python
'''converts genbank records into bed file usable for annotation''' 
import argparse
import sys;

from Bio.SeqFeature import CompoundLocation, FeatureLocation
from Bio import SeqIO
from pybedtools import BedTool, Interval;

from nrlbio.itertools_extension import flatten
from nrlbio.pybedtools_extension import construct_gff_interval


choice_types = ['misc_feature', 'STS', 'exon', 'mRNA', 'CDS', 'misc_RNA', 'gene']
default_types = ['mRNA', 'CDS', 'misc_RNA']


choice_miscrna_types = ['TEC', 'snRNA', 'lincRNA', 'unprocessed_pseudogene', 'antisense', 'retained_intron', 'transcribed_unprocessed_pseudogene', 'sense_intronic', 'processed_transcript', 'rRNA', 'transcribed_processed_pseudogene', 'miRNA', 'misc_RNA', 'snoRNA', 'sense_overlapping', 'processed_pseudogene'];
default_miscrna_types = ['TEC', 'snRNA', 'lincRNA', 'unprocessed_pseudogene', 'antisense', 'transcribed_unprocessed_pseudogene', 'sense_intronic', 'processed_transcript', 'rRNA', 'transcribed_processed_pseudogene', 'miRNA', 'misc_RNA', 'snoRNA', 'sense_overlapping', 'processed_pseudogene'];


parser = argparse.ArgumentParser(description='converts genbank records into bed file usable for annotation');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to genbank records");
parser.add_argument('-s', '--source', nargs = '?', default = "custom", help = "the name of the annotation source");
parser.add_argument('-t', '--rna_types', nargs = '+', default = default_types, choices = choice_types, type = str, help = "output only provided types of RNA. Please change default in case you really need it");
parser.add_argument('-mt', '--misc_rna_types', nargs = '+', default = default_miscrna_types, choices = choice_miscrna_types, type = str, help = "output only provided types of miscellaneous RNA");
args = parser.parse_args();



def get_exons(exons, start, end, strand):
		
	left_exons = [];
	right_exons = [];
	
	for e in exons:
		if(e.end <= start):
			left_exons.append(e);
		elif(e.start < start):
			left_exons.append(FeatureLocation(e.start, start, strand=strand))
		
		if(e.start >= end):
			right_exons.append(e);
		elif(e.end > end):
			right_exons.append(FeatureLocation(end, e.end, strand=strand));
	
	return left_exons, right_exons;
	

def get_introns(exons, start = None, end = None):
	if(not exons):
		return []
		
	strand = exons[0].strand
	fl = flatten([(x.start, x.end) for x in exons])
	introns = [ FeatureLocation(x[0], x[1], strand) for x in  zip(fl[1::2], fl[2::2]) ]
	
	if(start and start > exons[-1].end):
		introns.append(FeatureLocation(exons[-1].end, start, strand))
	if(end and end < exons[0].start):
		introns.insert(0, FeatureLocation(end, exons[0].start, strand))
		
	return introns	



def split_mrna(mrna, cds):
	strand = mrna.location.strand
	
	if(type(mrna.location) == CompoundLocation):
		exons = mrna.location.parts;
	else:
		exons = [mrna.location];
		
	if(type(cds.location) == CompoundLocation):
		cds_exons = cds.location.parts;
	else:
		cds_exons = [cds.location];		
		
	left_exons, right_exons = get_exons(exons, cds.location.start, cds.location.end, strand)

	left_introns = get_introns(left_exons, start = cds.location.start)
	right_introns = get_introns(right_exons, end = cds.location.end)
	cds_introns = get_introns(cds_exons)
				
	# we always return exons in an order 5'UTR, CDS, 3'UTR;		
	if(strand == 1):
		return left_exons, left_introns, cds_exons, cds_introns, right_exons, right_introns
	elif(strand == -1):
		return right_exons, right_introns, cds_exons, cds_introns, left_exons, left_introns
		
		
def location2gff(location, chrom, transcript_type, source, gene_id, regulation, transcription):
	if(location.strand == 1):
		strand = '+'
	elif(location.strand == -1):
		strand = '-'
	else:
		strand = '.'
	return construct_gff_interval(chrom, location.start+1, location.end, transcript_type, strand=strand, source=source, frame='.', attrs={"gene_id": gene_id, 'regulation': regulation, 'transcription': transcription})
		
		
def annotate_mrna(mrna, cds, chrom, source):
	utr5_exons, utr5_introns, cds_exons, cds_introns, utr3_exons, utr3_introns = split_mrna(mrna, cds);

	for location in cds_exons:
		print location2gff(location, chrom, 'mRNA', source, mrna.qualifiers['gene'][0], 'cds', 'exon');
	
	for location in cds_introns:
		print location2gff(location, chrom, 'mRNA', source, mrna.qualifiers['gene'][0], 'cds', 'intron');
		
	for location in utr5_exons:
		print location2gff(location, chrom, 'mRNA', source, mrna.qualifiers['gene'][0], 'utr5', 'exon');
	
	for location in utr5_introns:
		print location2gff(location, chrom, 'mRNA', source, mrna.qualifiers['gene'][0], 'utr5', 'intron');
		
	for location in utr3_exons:
		print location2gff(location, chrom, 'mRNA', source, mrna.qualifiers['gene'][0], 'utr3', 'exon');
	
	for location in utr3_introns:
		print location2gff(location, chrom, 'mRNA', source, mrna.qualifiers['gene'][0], 'utr3', 'intron');
		
		
def annotate_noncoding(feature, chrom, source):
	
	if(type(feature.location) == CompoundLocation):
		exons = feature.location.parts;
	else:
		exons = [feature.location];	
		
	introns = get_introns(exons);
	
	for location in exons:
		print location2gff(location, chrom, feature.type, source, feature.qualifiers['gene'][0], feature.qualifiers['note'][0], 'exon');
	
	for location in introns:
		print location2gff(location, chrom, feature.type, source, feature.qualifiers['gene'][0], feature.qualifiers['note'][0], 'intron');
		
	
	
def fetch_annotation(seq_record):
	for feature in seq_record.features:
		
		if(feature.type == "mRNA"):
			mrna = feature;
		elif(feature.type == "CDS"):
			annotate_mrna(mrna, feature, seq_record.name, args.source)
			
		elif(feature.type in args.rna_types):
			if(feature.type == "misc_RNA"):
				if(feature.qualifiers["note"][0] in args.misc_rna_types):
					annotate_noncoding(feature, seq_record.name, args.source)
				else:	
					pass;
			else:
				annotate_noncoding(feature, seq_record.name, args.source)



for path in args.path:
	for seq_record in SeqIO.parse(path, "genbank"):
		fetch_annotation(seq_record);



		
		
		