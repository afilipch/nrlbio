#! /usr/lib/python
'''Converts gff3 gene model into internal nrlbio genome system structure stored as yml file''' 
import argparse
import sys;
import yaml

from pybedtools import BedTool

from nrlbio.genome_system import gff3_to_genes






parser = argparse.ArgumentParser(description='Converts gff3 gene model into internal nrlbio genome system structure stored as yml file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the gene model, gff3 format");
parser.add_argument('--output', nargs = '?', required = True, type = str, help = "Path to the output");
args = parser.parse_args();

genes = [x.to_yaml() for x in gff3_to_genes(BedTool(args.path))]
with open(args.output, 'w') as f:
	f.write(yaml.dump(genes, default_flow_style=False))
	#for g in genes:
		#print g;
