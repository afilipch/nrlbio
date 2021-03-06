# /usr/bin/python
'''collections of classes and functions to deal with ensembl/ncbi/genbank records'''

import os, sys, re;
#import Bio;
#import pysam;
from collections import defaultdict;

from Bio.SeqFeature import CompoundLocation, FeatureLocation


def flanking_exons(exons, start, end, strand):
	'''selects for intervals(exons) provided only those which are not inside start and end
	
		exons list: list of Bio.SeqFeature.FeatureLocation instances to choose flanking ones from
		start int: all intervals on the left side of start will be selected
		end int: all intervals on the right side of end will be selected
		strand int: 1 for the top strand, -1 for the bottom strand, 0 if the strand is important but is unknown, or None if it does not matter.
	
	Return list: flanking intervals(Bio.SeqFeature.FeatureLocation instances);
	'''
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
	

adjusted_names = [str(x) for x in range(50)] + ['X', 'Y', 'MT'];
def adjust_name(name):
	'''adjust chromosome name from record to be consistent with chromosome names in other file formats'''
	if (name in adjusted_names):
		return "chr%s" % name
	else:
		return name;
		
	
	


def feature2fasta(feature, seq_record):
	'''parse Bio.SeqFeature into 2-line string representing single fasta entry
	
	feature Bio.SeqFeature: genbank feature to be converted into fasta entry
	seq_record Bio.SeqRecord: seq_record the feature refers to
	
	Return string: 2-line string representing single fasta entry
	'''
	try:
		key = feature.qualifiers['note'][-1];
	except:
		raise Exception('can not assign id to the feature%s\n' % feature)
	return ">%s\n%s" % (key, feature.extract(seq_record).seq);
	
	
def seq_record2fasta(seq_record, feature_types = ['mRNA', 'misc_RNA'], miscRNA_type = None):
	'''function converts every feature on Bio.SeqRecord into fasta entry and prints fasta file to stdout
	
	seq_record Bio.SeqRecord: seq_record consisting features to be converted into fasta entries
	feature_types list: function ouputs fasta entries only for those features whose Bio.SeqFeature.type is in the feature_types. If bool(feature_types)==False do not filter features according to Bio.SeqFeature.type. Bio.SeqFeature.type can be mRNA, misc_RNA, exon and etc
	miscRNA_type list: function ouputs fasta entries only for those features whose Bio.SeqFeature.qualifiers['note'] is in the miscRNA_type. If bool(miscRNA_type)==False do not filter features according to Bio.SeqFeature.qualifiers['note']. Bio.SeqFeature.qualifiers['note'] can be ncRNA, miRNA, snRNA ans etc.
	
	Return None: prints fasta file to stdout
	'''
	
	for feature in seq_record.features:
		if(not feature_types or feature.type in feature_types):
			if(feature.type == 'misc_RNA'):
				if (not miscRNA_type or feature.qualifiers['note'][0] in miscRNA_type):
					print(feature2fasta(feature, seq_record));
				else:
					pass;
			else:
				print(feature2fasta(feature, seq_record));
		else:
			pass;
	return None;
	
#def bed2feature()
	#my_start_pos = SeqFeature.ExactPosition(2)
	#my_end_pos = SeqFeature.ExactPosition(6)
	