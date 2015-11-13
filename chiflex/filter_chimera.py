#! /usr/lib/python
'''Filters chimeras on feature basis'''
import argparse
import sys;
import os;

import pysam;

from nrlbio.generators import generator_doublebed
from nrlbio.LRGFDR import lrg

import logging
# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
fh = logging.FileHandler(os.path.join("log", "chimera_filtering.txt"))
logger.addHandler(sh)
logger.addHandler(fh)


#get configuration
from nrlbio.config.config import load_config;
CONFIGURATION = load_config('lrg')



parser = argparse.ArgumentParser(description='Filters chimeras on feature basis');
parser.add_argument('-s', '--signal', nargs = '?', required = True, type = str, help = "path to the chimeras to be filtered, double gff format");
parser.add_argument('-c', '--control', nargs = '?', required = True, type = str, help = "path to the chimeras originated from decoy, double gff format");
parser.add_argument('-f', '--features', nargs = '+', default = ['AS1', 'AS2'], choices=['AS1', 'AS2', 'gap', 'minqstart'],  type = str, help = "Features to be used for filtering");
parser.add_argument('--fdr', nargs = '?', default = 0.05, type = float, help = "False Discovery Rate allowed");
#parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
args = parser.parse_args();


def as1(i1, i2):
	return int(i1.score);

def as2(i1, i2):
	return int(i2.score);

def gap(i1, i2):
	return int(i1.attrs['gap']);

def minqstart(i1, i2):
	return min(int(i1.attrs['qstart']), int(i2.attrs['qstart']));

features2indices = {'AS1': as1, 'AS2': as2, 'gap': gap, 'minqstart': minqstart}

def doublebed2list(i1, i2):
	return [features2indices[x](i1, i2) for x in args.features]

def list_generator(path):
	for i1, i2 in generator_doublebed(path):
		yield doublebed2list(i1, i2);
		
def apply_filter(path, filter_):
	for i1, i2 in generator_doublebed(path):
		x = doublebed2list(i1, i2);
		if(eval(filter_)):
			yield i1, i2
		else:
			pass; 
	



signal = list_generator(args.signal);
control = list_generator(args.control);
lrg_filter, rule, log_message = lrg(signal, control, entry='list', attributes = list(range(len(args.features))), attribute_names=args.features, fdr=args.fdr, **CONFIGURATION);

if(rule):
	for i1, i2 in apply_filter(args.signal, lrg_filter):
		sys.stdout.write(str(i1))
		sys.stdout.write(str(i2))	
	log_message = log_message
else:
	log_message = "Nothing passed the filtering\n"

logger.info(log_message)


