#! /usr/bin/python
'''Looks for the most abundant miRNA found in chimeras. Draws a histogram'''
import sys;
import argparse;
import os
from collections import defaultdict

from pybedtools import BedTool




parser = argparse.ArgumentParser(description='Looks for the most abundant miR7 targets. Draws a histogram');
parser.add_argument('path', metavar = 'N', nargs = '?', type = os.path.abspath, help = "Path to the interactions, single gff format(targets.gff)");
parser.add_argument('--norm', nargs = '?', default=False, const = True, type = bool, help = "If set, miRNAs expression values will be normalized");
parser.add_argument('--output', nargs = '?', default='mirexpr.svg', type = str, help = "Path to the output");
args = parser.parse_args();

mirids = defaultdict(float);

for interval in BedTool(args.path):
	mirids[interval.attrs['mirid']] += float(interval.attrs['n_uniq'])

top = list(sorted(mirids.items(), key = lambda x: x[1], reverse=True))[:10]
		
for k,v in top:
	sys.stderr.write("%s\t%d\n" % (k, v))
	
	
###PLOTING PART###	

values = [x[1] for x in top]
labels = [x[0] for x in top]

if(args.norm):
	nfactor = float(sum(mirids.values()))
	values = [x/nfactor for x in values]
	ylabel = "Fraction of chimeras"
else:
	ylabel = "Number of chimeras"

	
import matplotlib.pyplot as plt
#from nrlbio.pyplot_extension import remove_top_left_boundaries

plt.figure(1, figsize=(8,5))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
ax = plt.subplot(111)

plt.bar(range(len(values)), values, align='center', width=0.7, color='0.45', linewidth=0)
ax.set_xlim(-1, 10)
ax.set_xticks([x+0.2 for x in range(len(labels))])
ax.set_xticklabels( labels, rotation=45, ha='right')
ax.set_ylabel(ylabel)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
#ax.spines['left'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
for tic in ax.xaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
    
plt.savefig(args.output, bbox_inches='tight')  