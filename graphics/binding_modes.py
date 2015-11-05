# /usr/bin/python
'''Draws a plot of binding modes'''
import sys
import argparse
from collections import defaultdict
import math

from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np

from nrlbio.pyplot_extension import remove_top_left_boundaries
from nrlbio.mirna import fasta2mirnas, MODES_ORDER

parser = argparse.ArgumentParser(description='Draws a plot of binding modes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the interactiong RNA, double gff file with mode and mode_shuffled assigned");
parser.add_argument('--norm', nargs = '?', default=False, const = True, type = bool, help = "If set, control values for binding modes will renormalized, which removes a bias of binding modes order");
parser.add_argument('--output', nargs = '?', default='pattern.png', type = str, help = "Path to the output");
parser.add_argument('--title', nargs = '?', type = str, help = "Title for a plot");
args = parser.parse_args();

mode = defaultdict(int)
mode_shuffled = defaultdict(int)

bedtool = BedTool(args.path);
total = len(bedtool);

for interval in bedtool:
	mode[interval.attrs['mode']] += 1;
	mode_shuffled[interval.attrs['mode_shuffled']] += 1;
		
mode = np.array([mode[x] for x in MODES_ORDER])
mode_shuffled = np.array([mode_shuffled[x] for x in MODES_ORDER])

if(args.norm):
	coefs  = (float(total) + mode[0] - mode.cumsum())/(total + mode_shuffled[0] - mode_shuffled.cumsum())
	mode_shuffled = np.array([int(math.ceil(x[0]*x[1])) for x in zip(coefs, mode_shuffled)]);

	


#stop here

#plot section

#prepare an appropriate layout
plt.figure(1, figsize=(8,5))
ax = plt.subplot(111)
ax = remove_top_left_boundaries(ax)
plt.axis((0, len(pattern), 0, max(pattern)*1.2))
	
#plot real and control binding pattern
plt.plot(pattern, '0.2',  linewidth=2, label='interactions')
plt.plot(pattern_shuffled, '0.8',  linewidth=2, label='shuffled sequences')

#set labels and title
if(args.norm):
	plt.ylabel('base pairing [%]')
else:
	plt.ylabel('bases paired')	
plt.xlabel('miRNA position [nt]')
if(args.title):
	plt.title(args.title)
	
#set legend	
plt.legend(loc='upper right',prop={'size':10})

#output plot in PNG format
plt.savefig(args.output, bbox_inches='tight')