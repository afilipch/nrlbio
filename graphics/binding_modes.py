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
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the interacting RNA, double gff file with mode and mode_shuffled assigned");
parser.add_argument('--norm', nargs = '?', default=False, const = True, type = bool, help = "If set, binding modes will be normalized");
#parser.add_argument('--norm', nargs = '?', default=False, const = True, type = bool, help = "If set, control values for binding modes will renormalized, which removes a bias of binding modes order");
parser.add_argument('--output', nargs = '?', default='pattern.png', type = str, help = "Path to the output");
parser.add_argument('--title', nargs = '?', type = str, help = "Title for a plot");
args = parser.parse_args();

mode = defaultdict(int)
mode_shuffled = defaultdict(int)

bedtool = BedTool(args.path);
total = float(len(bedtool));

for interval in bedtool:
	mode[interval.attrs['mode']] += 1;
	mode_shuffled[interval.attrs['mode_shuffled']] += 1;
		
mode = np.array([mode[x] for x in MODES_ORDER])
mode_shuffled = np.array([mode_shuffled[x] for x in MODES_ORDER])

coefs  = (total + mode[0] - mode.cumsum())/(total + mode_shuffled[0] - mode_shuffled.cumsum())
mode_shuffled = np.array([int(math.ceil(x[0]*x[1])) for x in zip(coefs, mode_shuffled)]);

if(args.norm):
	mode = mode/total*100
	mode_shuffled = mode_shuffled/total*100
else:
	mode = mode/2;
	mode_shuffled = mode_shuffled/2
	
#print mode
#print mode_shuffled
#sys.exit();
	


#stop here

#plot section

#prepare an appropriate layout
plt.figure(1, figsize=(8,5))
ax = plt.subplot(111)
remove_top_left_boundaries(ax)
plt.axis((0, len(mode)+1, 0, max(mode)*1.2))

#set bins and boundaries
boundaries = range(0, len(mode));
bins = range(0, len(mode)+1);

#plot real and control binding pattern
plt.hist((boundaries,boundaries), weights=(mode,mode_shuffled), bins=bins, label=('interactions(%d)' % (total/2), 'shuffled sequences'), align='right', rwidth=0.7, color=('0.2', '0.8'))


#set labels and title
plt.xlabel('binding mode')	
if(args.norm):
	plt.ylabel('fraction of interactions [%]')
else:
	plt.ylabel('number of interactions')	
if(args.title):
	plt.title(args.title)
	
#set xlabels	
plt.xticks(range(1, len(mode)+1));
ax.set_xticklabels(MODES_ORDER, rotation=0)
	
#set legend	
plt.legend(loc='upper right',prop={'size':10})


#output plot in PNG format
plt.savefig(args.output, bbox_inches='tight')