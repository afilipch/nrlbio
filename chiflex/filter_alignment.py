#! /usr/lib/python
'''filters mapping results on read features basis'''
import argparse
import os;

import pysam;

from nrlbio import sam_statistics
from nrlbio.samlib import filter_generator, apply_filter
from nrlbio.LRGFDR import lrg

parser = argparse.ArgumentParser(description='filters mapping results on read features basis');
parser.add_argument('-s', '--signal', nargs = '?', required = True, type = str, help = "path to the sam file to be filtered");
parser.add_argument('-c', '--control', nargs = '?', required = True, type = str, help = "path to the sam file originated from decoy");
parser.add_argument('-f', '--features', nargs = '+', default = ['AS'], type = str, help = "read features to be used for filtering");
parser.add_argument('-o', '--output', nargs = '?', default = "sam", type = str, help = "path to the output folder");
parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for output files, should reflect nature of mapping reference");
args = parser.parse_args();



signal = filter_generator(pysam.Samfile(args.signal), args.features);
control = filter_generator(pysam.Samfile(args.control), args.features);
lrg_filter, rule = lrg(signal, control, entry='list', attribute_names=args.features, support = 0.02, maxiter = 20,  fdr=0.1, lookforward=10, ncsupport=0.1, nciter=1)

samfile = pysam.Samfile(args.signal)
filtered = pysam.Samfile(os.path.join(args.output, "%s.filtered.bam" % args.name), "wb", template=samfile)
for ar in apply_filter(samfile, args.features, lrg_filter):
	filtered.write(ar)


