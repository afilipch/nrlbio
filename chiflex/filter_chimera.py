#! /usr/lib/python
'''filters mapping results on read features basis'''
import argparse
import sys;

import pysam;

from nrlbio.chimera import filter_generator, filter_generator_doublebed, apply_filter, apply_filter_doublebed
from nrlbio.LRGFDR import lrg

parser = argparse.ArgumentParser(description='filters mapping results on read features basis');
parser.add_argument('-s', '--signal', nargs = '?', required = True, type = str, help = "path to the sam file to be filtered");
parser.add_argument('-c', '--control', nargs = '?', required = True, type = str, help = "path to the sam file originated from decoy");
parser.add_argument('-f', '--features', nargs = '+', default = ['AS1', 'AS2'], type = str, help = "read features to be used for filtering");
parser.add_argument('--fdr', nargs = '?', default = 0.05, type = float, help = "False Discovery Rate allowed");
#parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
#parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for filtered sam");
args = parser.parse_args();

features2indices = {'AS1': 4, 'AS2': 10, 'gap': 13, 'qstart1': 14, 'qstart2': 17, 'bias1': 16, 'bias2': 19}
indices = [];
for f in args.features:
	indices.append(features2indices[f]);



signal = filter_generator_doublebed(args.signal, indices);
control = filter_generator_doublebed(args.control, indices);
lrg_filter, rule = lrg(signal, control, entry='list', attributes = indices, attribute_names=args.features, support = 0.02, maxiter = 20,  fdr=args.fdr, lookforward=10, ncsupport=0.1, nciter=2)

if(not rule):
	samfile.close()
	filtered.close()
	sys.exit("Nothing passed the filtering\n")

for l in apply_filter_doublebed(args.signal, indices, lrg_filter):
	print l





