# /usr/bin/python
'''Draws plots for sensitivity and specifity of chiflex as a function of mismatch rate for different conditions (background model, length of reads)'''
import sys
import os
import argparse
from collections import defaultdict
import math


import matplotlib.pyplot as plt
import numpy as np

from nrlbio.pyplot_extension import remove_top_left_boundaries


parser = argparse.ArgumentParser(description='Draws plots for sensitivity and specifity of chiflex as a function of mismatch rate for different conditions (background model, length of reads)');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the sensitivity and specificity tables. Tab delimited format");
parser.add_argument('--names', nargs = '+', required=True,type = str, help = "Names of conditions. Each condition has to correspond to the input file");
parser.add_argument('--output', nargs = '?', required=True, type = str, help = "Path to the output directory for plots");
parser.add_argument('--title', nargs = '?', type = str, help = "Title for a plot");
parser.add_argument('--xlabel', nargs = '?', default = 'mrate', choices=('mrate', 'length'), type = str, help = "Type of xlabel");
args = parser.parse_args();

tables = defaultdict(dict)

for path, name in zip(args.path, args.names):
	with open(path) as f:
		mmrate = f.next().strip().split("\t")[1:]
		for l in f:
			a = l.strip().split("\t")
			tables[name][a[0]] = np.array([float(x) for x in a[1:]])



	
#print mmrate
#for k, v in tables['shuffled'].items():
	#print k, v
#sys.exit();
	

##############################################################################
#Plotting

def hist(data1, data2, output, labels, xticklabels, title, xlabel):
	#prepare an appropriate layout
	plt.figure(1, figsize=(8,5))
	ax = plt.subplot(111)
	remove_top_left_boundaries(ax)
	plt.axis((0, len(data1)+1, 0, max(max(data2), max(data1))*1.2))

	#set bins and boundaries
	boundaries = range(0, len(data1));
	bins = range(0, len(data2)+1);

	#plot real and control binding pattern
	plt.hist((boundaries,boundaries), weights=(data1, data2), bins=bins, label=labels, align='right', rwidth=0.7, color=('0.1', '0.5'))


	#set labels and title
	if(xlabel == 'mrate'):
		plt.xlabel('mismatch rate [%]')
	elif(xlabel == 'length'):
		plt.xlabel('read length [nt]')
	plt.ylabel('fraction of reads')
	plt.title(title)
		
	#set xlabels	
	plt.xticks(range(1, len(data1)+1));
	ax.set_xticklabels(xticklabels, rotation=0)
		
	#set legend	
	plt.legend(loc='upper right',prop={'size':10})


	#output plot in PNG format
	plt.savefig(output, bbox_inches='tight')
	plt.close();
	
	

for name, table in tables.items():
	#plot sensitivity
	sens_chimera_mapping = table['sensitivity_chimera_(mapping)']
	sens_chimera_total = table['total_sensitivity_chimera']	
	hist(sens_chimera_mapping, sens_chimera_total, os.path.join(args.output, "sens_chimera_%s.png" % name), labels=('mapping', 'filtering'), xticklabels=mmrate, title="sensitivity of chimera recovery", xlabel=args.xlabel)
	
	sens_single_mapping = table['sensitivity_single_(mapping)']
	sens_single_total = table['total_sensitivity_single']
	hist(sens_single_mapping, sens_single_total, os.path.join(args.output, "sens_single_%s.png" % name), labels=('mapping', 'filtering'), xticklabels=mmrate, title="sensitivity of single reads recovery", xlabel=args.xlabel)
	
	
	#plot specificity
	spec_chimera_mapping = table['specificity_chimera_(mapping)']
	spec_chimera_total = table['specificity_chimera']	
	hist(spec_chimera_mapping, spec_chimera_total, os.path.join(args.output, "spec_chimera_%s.png" % name), labels=('mapping', 'filtering'), xticklabels=mmrate, title="specificity of chimera recovery", xlabel=args.xlabel)
	
	spec_single_mapping = table['specificity_single_(mapping)']
	spec_single_total = table['specificity_single']
	hist(spec_single_mapping, spec_single_total, os.path.join(args.output, "spec_single_%s.png" % name), labels=('mapping', 'filtering'), xticklabels=mmrate, title="specificity of single reads recovery", xlabel=args.xlabel)
	
	
	#plot fdr
	fdr_single_real = 100 - table['specificity_single']
	fdr_single_estimated = table['estimated_fdr_single']
	hist(fdr_single_real, fdr_single_estimated, os.path.join(args.output, "fdr_single_%s.png" % name), labels=('real', 'estimated'), xticklabels=mmrate, title="FDR of single reads recovery", xlabel=args.xlabel)
	
	fdr_chimera_real = 100 - table['specificity_chimera']
	fdr_chimera_estimated = table['estimated_fdr_chimera']
	hist(fdr_chimera_real, fdr_chimera_estimated, os.path.join(args.output, "fdr_chimera_%s.png" % name), labels=('real', 'estimated'), xticklabels=mmrate, title="FDR of chimera recovery", xlabel=args.xlabel)
	
	
	
	
	
	
	
	