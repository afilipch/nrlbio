# /usr/bin/python
'''support genomic annotion of bed and doublebed files'''
import sys;
import os;
from collections import defaultdict, OrderedDict;

from nrlbio.numerictools import maxes 

STAT_ATTRS = 'transcript_types', 'transcription_blocks', 'mRNA', 'regulation_regions', 'misc_RNA'

#order of features to select fo report
order = {}
order['type'] = {'mRNA':1, 'misc_RNA':2, 'regulation':2, 'intergenic':3};

order['regulation'] = {};
order['regulation']['mRNA'] = {'utr5':1, 'utr3':1, 'cds':1}
order['regulation']['regulation'] = {'TF_binding_site':1, 'enhancer':2, 'promoter':2, 'promoter_flanking_region': 3, 'open_chromatin_region':4, 'CTCF_binding_site':4}

order['regulation']['misc_RNA'] = {"lincRNA":1, "snRNA":1, "miRNA":1, 'known_ncrna':1, "pseudogene":1, 'TEC':1, 
'transcribed_unprocessed_pseudogene':2, 'unprocessed_pseudogene':2, 'processed_pseudogene':2, "transcribed_processed_pseudogene":2,
'non_coding':3, 'retained_intron':3,
'antisense':4, 'sense_intronic':4, 'sense_overlapping':4, '3prime_overlapping_ncrna':4,
'misc_RNA':5,
'processed_transcript':6}

order["transcription"]={'exon':1, 'intron':1}


def interval2annotation(interval):
	ann = defaultdict(lambda: defaultdict(list));
	for i, t in enumerate(interval.attrs['type'].split(",")):
		for j, r in enumerate(interval.attrs['regulation'].split(",")[i].split("|")):
			ann[t][r] = interval.attrs['transcription'].split(",")[i].split("|")[j].split('%');
	return ann;
	
	
def assign_stat_attributes(interval):
	transcription_blocks = set()
	types_reg = {};
	ann = interval2annotation(interval);
	transcript_types, stub = maxes(ann.keys(), key_function=lambda x: -1*order["type"][x])
	for t in transcript_types:
		if(any(ann[t].keys())):
			types_reg[t], stub = maxes(ann[t].keys(), key_function=lambda x: -1*order["regulation"][t][x]);
			for r in types_reg[t]:
				transcription_blocks.update(ann[t][r]);
				
	interval.attrs['mRNA'] = ",".join(types_reg.get('mRNA',[]))
	interval.attrs['regulation_regions'] = ",".join(types_reg.get('regulation', []))
	interval.attrs['misc_RNA'] = ",".join(types_reg.get('misc_RNA', []))
	
	interval.attrs['transcript_types'] = ",".join(transcript_types);
	interval.attrs['transcription_blocks'] = ",".join(transcription_blocks);
	


def annotation_generator(bedtool, attributes):
	temp = defaultdict(set)
	for interval in bedtool:
		for c , t in enumerate(interval.attrs['type'].split(",")):
			temp[t].add(interval.attrs['regulation'].split(",")[c])
		yield interval;
	for k,v in temp.items():
		print k, v


