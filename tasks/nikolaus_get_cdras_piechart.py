#! /usr/bin/python
'''Draws a circular CDR1as plot with discovered chimeric targets''' 
import argparse
import sys;
from collections import defaultdict
from itertools import chain
import math

import numpy as np
import matplotlib.pyplot as plt

from nrlbio.generators import generator_doublebed



parser = argparse.ArgumentParser(description='Draws a circular CDR1as plot with discovered chimeric targets');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the mirna:target interactions to select from, double gff");
parser.add_argument('--mirids', nargs = 2, required = True, type = str, help = "miRNA id to select");
parser.add_argument('--cid', nargs = '?', required = True, type = str, help = "cicular id to select");
parser.add_argument('--length', nargs = '?', required = True, type = float, help = "length of circle");

args = parser.parse_args();

mirids = args.mirids
cid = args.cid
dintervals = defaultdict(list)

for i1, i2 in generator_doublebed(args.path):
	if(i1.chrom in mirids and i2.chrom == cid):
		dintervals[i1.chrom].append(i2)



#set up the plot

bottom = 8
max_height = 8
pos2grade = (2*np.pi)/args.length

ax = plt.subplot(111, polar=True)
ax.grid(False)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.spines['polar'].set_visible(False)



#fill the plot

colors = ['skyblue', 'coral']
maxcov = max([float(x.attrs['n_uniq'])/len(x) for x in  chain.from_iterable(dintervals.values())])


lintervals = [dintervals[x] for x in mirids]

for intervals, color in zip(lintervals, colors):
	intervals.sort(key= lambda x: x.start)
	#for i in intervals:
		#print i.start, i.end
	#print	
	theta = [-(x.end*pos2grade)+np.pi/2 for x in intervals]
	width = [len(x)*pos2grade for x in intervals]
	coverage = [float(x.attrs['n_uniq'])/len(x) for x in intervals]
	if(color=='skyblue'):
		radii = [x*max_height/maxcov for x in coverage]
		bars = ax.bar(theta, radii, width=width, bottom=bottom, color=color, linewidth=0)
	else:
		radii = [x*max_height/maxcov for x in coverage]
		bars = ax.bar(theta, radii, width=width, bottom=np.repeat(bottom,len(intervals))-radii, color=color, linewidth=0)
		
		
legend = plt.legend(['miR-7', 'miR-671'], loc=(0.9, 0.9), frameon=False)
ax.plot(np.linspace(0, 2*np.pi, 100), [bottom]*100, color='r', linewidth=3)


plt.show()











