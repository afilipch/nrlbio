#! /usr/bin/python
'''Draws a circular CDR1as plot with discovered chimeric targets and conservation''' 
import argparse
import sys;
#from collections import defaultdict
#from itertools import chain
from math import log

from pybedtools import BedTool
import numpy as np
import matplotlib.pyplot as plt





parser = argparse.ArgumentParser(description='Draws a circular CDR1as plot with discovered chimeric targets and conservation');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "path to the coverage for miRNAs of interest (mir-7, mir-671), bedgraph format");
parser.add_argument('--output', nargs = '?', type = str, help = "Path to the output");
parser.add_argument('--cid', nargs = '?', required = True, type = str, help = "cicular id to select");
parser.add_argument('--length', nargs = '?', required = True, type = float, help = "length of circle");
args = parser.parse_args();

def get_coverage(intervals):
	return [log(float(x[3])+1, 2) for x in intervals]


lintervals = [ [], [] ]

for interval in BedTool(args.path[0]):
	if(interval.chrom == args.cid):
		lintervals[0].append(interval);
	
for interval in BedTool(args.path[1]):
	if(interval.chrom == args.cid):
		lintervals[1].append(interval);
		
		
		
#set up the plot
bottom = 6
max_height = 5
pos2grade = (2*np.pi)/args.length

ax = plt.subplot(111, polar=True)
ax.grid(False)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.spines['polar'].set_visible(False)


#fill the plot
colors = ['skyblue', 'coral'];
maxcov = max(get_coverage(lintervals[0]) + get_coverage(lintervals[1]))


for intervals, color in zip(lintervals, colors):
	intervals.sort(key= lambda x: x.start)

	theta = [-(x.end*pos2grade)+np.pi/2 for x in intervals]
	width = [len(x)*pos2grade for x in intervals]
	coverage = get_coverage(intervals)
	
	radii = [x*max_height/maxcov for x in coverage]
	#bars = ax.bar(theta, radii, width=width, bottom=bottom, color=color, linewidth=0)
	
	if(color=='skyblue'):
		bars = ax.bar(theta, radii, width=width, bottom=bottom, color=color, linewidth=0)
	else:
		bars = ax.bar(theta, radii, width=width, bottom=np.repeat(bottom,len(intervals))-radii, color=color, linewidth=0)
		
		
legend = plt.legend(['miR-7', 'miR-671'], loc=(0.9, 0.9), frameon=False)
ax.annotate("splice junction", xy=(np.pi/2, bottom), xytext=(np.pi/2, bottom+3), arrowprops=dict(arrowstyle="->"), ha='center')
ax.plot(np.linspace(0, 2*np.pi, 100), [bottom]*100, color='r', linewidth=1)

if(args.output):
	plt.savefig(args.output)
else:
	plt.show()
	





