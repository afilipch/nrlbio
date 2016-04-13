#! /usr/bin/python
'''Draws a circular CDR1as plot with discovered chimeric targets and conservation''' 
import argparse
import sys;
#from collections import defaultdict
#from itertools import chain
#import math

from pybedtools import BedTool
import numpy as np
import matplotlib.pyplot as plt





parser = argparse.ArgumentParser(description='Draws a circular CDR1as plot with discovered chimeric targets and conservation');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the coverage for miRNA of interest, bedgraph format");
parser.add_argument('--conservation', nargs = '?', required = True, type = str, help = "path to the conservation score, bedgraph format");
parser.add_argument('--cid', nargs = '?', required = True, type = str, help = "cicular id to select");
parser.add_argument('--length', nargs = '?', required = True, type = float, help = "length of circle");
args = parser.parse_args();


lintervals = [ [], [] ]

for interval in BedTool(args.path):
	if(interval.chrom == args.cid):
		lintervals[0].append(interval);
	
for interval in BedTool(args.conservation):
	if(interval.chrom == args.cid):
		lintervals[1].append(interval);
		
		
		
#set up the plot
bottom = 8
max_height = 4
pos2grade = (2*np.pi)/args.length

ax = plt.subplot(111, polar=True)
ax.grid(False)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.spines['polar'].set_visible(False)


#fill the plot
colors = ['skyblue', 'coral'];

for intervals, color in zip(lintervals, colors):
	print 'n'
	intervals.sort(key= lambda x: x.start)

	theta = [-(x.end*pos2grade)+np.pi/2 for x in intervals]
	width = [len(x)*pos2grade for x in intervals]
	coverage = [max(float(x[3]), 0) for x in intervals]
	maxcov = max(coverage)
	
	radii = [x*max_height/maxcov for x in coverage]
	
	if(color=='skyblue'):
		bars = ax.bar(theta, radii, width=width, bottom=bottom, color=color, linewidth=0)
	else:
		bars = ax.bar(theta, radii, width=width, bottom=np.repeat(bottom,len(intervals))-radii, color=color, linewidth=0)
		
		
legend = plt.legend(['coverage', 'conservation'], loc=(0.9, 0.9), frameon=False)
ax.annotate("splice junction", xy=(np.pi/2, 8), xytext=(np.pi/2, 11), arrowprops=dict(arrowstyle="->"), ha='center')
#plt.arrow(0, 8, 0, -8)
ax.plot(np.linspace(0, 2*np.pi, 100), [bottom]*100, color='r', linewidth=3)


plt.show()





