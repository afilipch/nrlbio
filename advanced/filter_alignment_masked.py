#! /usr/lib/python
'''filters mapping results on read features basis'''
import argparse
import os;
import sys;

import pysam;

from nrlbio.samlib import filter_generator, apply_filter, get_attributes_masked
from nrlbio.LRGFDR import lrg

import logging
# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
fh = logging.FileHandler(os.path.join("log", "alignment_filtering.txt"))
logger.addHandler(sh)
logger.addHandler(fh)

#get configuration
from nrlbio.config.config import load_config;
CONFIGURATION = load_config('lrg')

parser = argparse.ArgumentParser(description='Filters mapping results on read features basis');
parser.add_argument('-s', '--signal', nargs = '?', required = True, type = str, help = "path to the sam file to be filtered");
parser.add_argument('-c', '--control', nargs = '?', required = True, type = str, help = "path to the sam file originated from decoy");
parser.add_argument('-f', '--features', nargs = '+', default = ['AS'], type = str, help = "read features to be used for filtering");
parser.add_argument('--fdr', nargs = '?', default = 0.05, type = float, help = "False Discovery Rate allowed");
parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
parser.add_argument('-n', '--output', nargs = '?', required = True, type = str, help = "Path for the filtered sam/bam file");
args = parser.parse_args();


ss = pysam.Samfile(args.signal)
sc = pysam.Samfile(args.control)
filtered = pysam.Samfile(args.output, "wb", template=ss)


try: 
	ss.next()
	ss.reset()
except:	
	sys.exit("Error: bam/sam signal file is empty\n")
try: 
	sc.next()
	sc.reset()
except:	
	sys.stderr.write("Warning: bam/sam signal file is empty\nAll entries of signal file will pass the filter\n");
	for ar in ss:
		filtered.write(ar)
	sys.exit();	
	

signal = filter_generator(ss, args.features, ga=get_attributes_masked);
control = filter_generator(sc, args.features, ga=get_attributes_masked);
lrg_filter, rule, log_message = lrg(signal, control, entry='list', attribute_names=args.features, fdr=args.fdr, **CONFIGURATION)

ss.close();
sc.close();



samfile = pysam.Samfile(args.signal);


if(not rule):
	samfile.close()
	filtered.close()
	sys.exit("Nothing passed the filtering\n")

for ar in apply_filter(samfile, args.features, lrg_filter, ga=get_attributes_masked):
	filtered.write(ar)
samfile.close()	



logger.info(log_message)








		

	