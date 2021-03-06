# /usr/bin/python
'''Draws a plot of nontemplate addiction to miRNAs for wt versus ko'''
import sys
import argparse
from collections import defaultdict
import math

from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np

from nrlbio.pyplot_extension import remove_top_left_boundaries


parser = argparse.ArgumentParser(description='Draws a plot of nontemplate addiction to miRNAs for wt versus ko');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "Path to the output of 'nontemplate_advanced.py' script, WT and KO");
parser.add_argument('--norm', nargs = '?', default=False, const = True, type = bool, help = "If set, binding modes will be normalized");
parser.add_argument('--output', nargs = '?', type = str, help = "Path to the output");
parser.add_argument('--extlen', nargs = '?', required=True, type = int, help = "The length of downstream extensions oof mature miRNAs");
args = parser.parse_args();

colors = ('greenyellow', 'skyblue', 'peachpuff', 'navy')
titles = ('WT', 'KO')

def path2weights(path):
	counts = [[], [], [], []]
	xlabels = [];
	totals = [];
	
	with open(path) as f:
		f.next();
		for l in f:
			a = l.strip().split('\t')
			for i in range(4):
				counts[i].append(float(a[i+1]))
				
			xlabels.append("\n".join(a[5].split(',')));
			totals.append(sum([float(a[x+1]) for x in range(4)]))
			
	maximum = max(totals);
	length = len(counts[1])

	xlabels.insert(length-args.extlen, ' ')
	for l in counts:
		l.insert(length-args.extlen, 0)
	length+=1;	

	counts = [np.array(x)/maximum for x in counts]

	return counts, xlabels, length


def draw_subplot(ax, counts, xlabels, length, title):
	remove_top_left_boundaries(ax)
	ylim = 1.01
	ax.set_ylim((0,ylim))
	ax.set_title(title)
	#plt.axis((0, length+1, 0, ylim))

	#set bins and boundaries
	boundaries = range(0, length);
	bins = range(0, length+1);

	#plot real and control binding pattern
	ax.hist([boundaries]*4, weights=counts, bins=bins, label=('A', 'C', 'T', 'G'), align='right', rwidth=0.7, color=colors, stacked=True)

	ax.set_xticks(range(1, length+1));
	ax.set_xticklabels(xlabels, rotation=0)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	for tic in ax.xaxis.get_major_ticks():
		tic.tick1On = tic.tick2On = False

	ax.set_ylabel('fraction of mappings')

	ax.legend(loc=(0.85, 0.75),prop={'size':10}, frameon=False)	


fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False, figsize=(16, 5))
for ax, path, title in zip(axes, args.path, titles):
	counts, xlabels, length = path2weights(path);
	draw_subplot(ax, counts, xlabels, length, title);


#output plot
if(args.output):
	plt.savefig(args.output, bbox_inches='tight')
else: 
	plt.show();