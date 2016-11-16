#! /usr/bin/python
'''Extracts pseudogenes from ensembl annotation'''


import sys;
import argparse

from pybedtools import BedTool

from nrlbio.pybedtools_extension import gff2bed


parser = argparse.ArgumentParser(description='Extracts pseudogenes from ensembl annotation');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome system annotation, gff3 ensembl format");
args = parser.parse_args();


exons = [];
btypes = set(['unitary_pseudogene', 'IG_C_pseudogene', 'polymorphic_pseudogene', 'TR_J_pseudogene', 'IG_J_pseudogene', 'TEC', 'TR_V_pseudogene', 'IG_V_pseudogene', 'pseudogene', 'unprocessed_pseudogene', 'transcribed_unprocessed_pseudogene', 'translated_unprocessed_pseudogene', 'processed_transcript', 'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene', 'processed_pseudogene'])

for interval in BedTool(args.path):
	if(interval.attrs.get('ID', '').split(':')[0] == 'gene'):
		if(interval.attrs['biotype'] in btypes):
			sys.stdout.write(str(gff2bed(interval)))
		
