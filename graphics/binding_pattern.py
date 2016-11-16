# /usr/bin/python
'''Draws a plot of binding pattern'''
import sys
import argparse
from collections import defaultdict

from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np

from nrlbio.pyplot_extension import remove_top_left_boundaries


parser = argparse.ArgumentParser(description='Draws a plot of binding pattern');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the interactiong RNA, double gff file with pattern and pattern_shuffled assigned");
parser.add_argument('--norm', nargs = '?', default=False, const = True, type = bool, help = "If set, binding pattern will be normalized");
parser.add_argument('--length', nargs = '?', default=0, type = int, help = "Length of binding pattern, if not set binding pattern will be calculated for the longest sRNA in interactions");
parser.add_argument('--output', nargs = '?', type = str, help = "Path to the output");
parser.add_argument('--title', nargs = '?', type = str, help = "Title for a plot");
parser.add_argument('--control', nargs = '?', default = 'shuffled sequences', type = str, help = "Type of control: shuffled sequences or shuffled interactions");
args = parser.parse_args();

pattern = defaultdict(float)
pattern_shuffled = defaultdict(float)

bedtool = BedTool(args.path);
total = len(bedtool);

for interval in bedtool:
	for c,p in enumerate(interval.attrs['pattern'].split(",")):
		pattern[c] += float(p);
	for c,p in enumerate(interval.attrs['pattern_shuffled'].split(",")):
		pattern_shuffled[c] += float(p);
		
pattern = np.array([pattern[x] for x in sorted(pattern.keys())])
pattern_shuffled = np.array([pattern_shuffled[x] for x in sorted(pattern_shuffled.keys())])

if(args.norm):
	pattern = pattern/total*100;
	pattern_shuffled = pattern_shuffled/total*100;
else:
	pattern = pattern/2
	pattern_shuffled = pattern_shuffled/2
	
	
if(args.length):
	pattern = pattern[:args.length]
	pattern_shuffled = pattern_shuffled[:args.length]



#plot section
colors = ('0.2', '0.8')

#prepare an appropriate layout
plt.figure(1, figsize=(8,5))
ax = plt.subplot(111)
remove_top_left_boundaries(ax)
if(args.norm):
	plt.axis((0,len(pattern),0, 100))
else:
	plt.axis((0, len(pattern), 0, max(pattern)*1.2))
	
#plot real and control binding pattern
plt.plot(pattern, colors[0],  linewidth=2, label='interactions(%d)' % (total/2))
plt.plot(pattern_shuffled, colors[1],  linewidth=2, label=args.control)

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
if(args.output):
	plt.savefig(args.output, bbox_inches='tight')
else: 
	plt.show();

	
	
	


