#! /usr/lib/python
'''Generates fasta file from given markov models'''

import argparse
import os
import sys

from Bio import SeqIO

from nrlbio.HMM import MultiMarkov
from nrlbio.systemtools import onlyfiles

parser = argparse.ArgumentParser(description='Generates fasta file from given markov models');
parser.add_argument('--models', nargs = '?', default = '', type = str, help = "path to the folder with markov models. For each serialized markov model(yml file) found in the directory fasta record will be generated. Length of the record is derived from model length, name is derived from file name");
parser.add_argument('--global_model', nargs = '?', default = '', type = str, help = "path to the serialized markov model. This only one model will be used to generate many fasta entries defined via \'--names\' and \'--lengths\' attributes. Ignored if \'--models\' is set");
parser.add_argument('--names', nargs = '+', default=[], type = str, help = "Names of the generated sequence record. NOTE: \'random_\' prefix will be added automatically. The order of arguments should correspond to the one of  \'--lengths\' arguments. Ignored if \'--models\' is set");
parser.add_argument('--lengths', nargs = '+', default=[], type = str, help = "Lengths of the generated sequence records. The order of arguments should correspond to the one of  \'--names\' arguments. Ignored if \'--models\' is set");
parser.add_argument('-ll', '--line_length', nargs = '?', default = 100, type = int, help = "Max length of line in output fasta file. Makes it more readable");
args = parser.parse_args();


if(args.models):
	for fname in onlyfiles(args.models):
		header, ext = os.path.splitext(fname);
		serialized = os.path.join(args.models, fname)
		if(ext == ".yml"):
			mm = MultiMarkov.deserialize(serialized);
			try:
				mm = MultiMarkov.deserialize(serialized);
			except:
				sys.stderr.write('Warning! %s file cannot be converted to MultiMarkov object. It is skipped then\n\n' % serialized)
				mm = None;
				
			if(mm):
				print ">%s" % header
				k = 0
				for s in mm.generate_string(args.line_length):
					k += s;
				print k;	
			else:
				pass;
		else:
			pass;
	
elif(args.global_model):
	try:
		mm = MultiMarkov.deserialize(args.global_model);
	except:
		sys.exit('Warning! %s file cannot be converted to MultiMarkov object' % args.global_model);
		
	for header, length in zip(args.names, args.length):
		print ">%s" % header
		for s in mm.generate_string(args.line_length, length=length):
			print s;	
else:
	sys.exit('\'--models\' or \'--global_model\' has to be set\n\n')

