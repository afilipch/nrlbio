#! /usr/lib/python
'''filters mapping results on read features basis'''
import argparse
import os;

import pysam;

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
lrg_filter, rule, signal_total, control_total, support_total, fdr_total = lrg(signal, control, entry='list', attribute_names=args.features, support = 0.02, maxiter = 20,  fdr=0.1, lookforward=10, ncsupport=0.1, nciter=1)

print lrg_filter
print rule
print 
print signal_total, control_total, support_total, fdr_total

#signal_real, control_real, total_real = 0,0,0

samfile = pysam.Samfile(args.signal);
filtered = pysam.Samfile(os.path.join(args.output, "%s.filtered.bam" % args.name), "wb", template=samfile)
for ar in apply_filter(samfile, args.features, lrg_filter):
	#signal_real += 1;
	filtered.write(ar)
samfile.close()	



#control_samfile = pysam.Samfile(args.control);
#for ar in apply_filter(control_samfile, args.features, lrg_filter):
	#control_real += 1;
	
	
#samfile = pysam.Samfile(args.signal);
#for aligned_read in samfile.fetch(until_eof=True):
	#if(not aligned_read.is_unmapped):
		#total_real+=1
		
		
#support_real = signal_real/(total_real+0.1)
#fdr_real = control_real/(signal_real+control_real+0.1)
#print 
#print signal_real, control_real, support_real, fdr_real
		

	