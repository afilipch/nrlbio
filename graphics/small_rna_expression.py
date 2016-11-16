#! /usr/bin/python
'''Draws scatter plot of small RNA differential expression'''
import argparse;
import sys;
import os;
import matplotlib.pyplot as plt
import matplotlib.lines as plines

from math import log
from collections import defaultdict

parser = argparse.ArgumentParser(description='Draws scatter plot of small RNA differential expression');
parser.add_argument('path', metavar = 'N', nargs = '?', type = os.path.abspath, help = "Path to the mirna differential expression, tsv file)");
parser.add_argument('--sizes', nargs = 2, required = True, type = int, help = "Number of WT and KO samples");
parser.add_argument('--output', nargs = '?', default='', type = str, help = "Path to the graphical output");
args = parser.parse_args();

size1, size2 = args.sizes;
ylims = [-3.5,3.5]
selected_labels = {'mmu-miR-7a-5p': ('b', 'o'), 'mmu-miR-7b-5p': ('g', 'o'), 'mmu-miR-7a-3p': ('b', '*'), 'mmu-miR-7b-3p': ('g', '*'), 'mmu-miR-671-5p': ('y', 'o'), 'mmu-miR-671-3p': ('y', '*')}

texpr = defaultdict(lambda: defaultdict(list))
with open(args.path) as f:
	f.next();
	for line in f:
		l = line.strip().split("\t")
		expr1 = [float(x) for x in l[1:1+size1]]
		expr2 = [float(x) for x in l[1+size1:1+size1+size2]]
		tm = l[0].split('-')
		if(len(tm)>4):
			mirid = '-'.join(tm[:3] + tm[-1:])
		else:
			mirid = l[0]	
		texpr[mirid][0].extend(expr1)
		texpr[mirid][1].extend(expr2)
		
mirexpr = {};
for mirid, d in texpr.items():
	expr1 = sum(d[0])/len(d[0])
	expr2 = sum(d[1])/len(d[1])
	lfc = log(expr2/expr1, 2)
	mirexpr[mirid] = (log(expr1, 2), lfc)



plt.figure(1, figsize=(8,5))
#plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
ax = plt.subplot(111)
plt.ylim(*ylims)
plt.xlim(0, 20)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
#ax.spines['left'].set_visible(False)

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
for tic in ax.xaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
for tic in ax.yaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False    


plt.plot([0,20], [0,0], 'k')
plt.plot([0,20], [0.5,0.5], 'k--')
plt.plot([0,20], [-0.5,-0.5], 'k--')

selected_plines = []

for label, (expr, lfc) in mirexpr.iteritems():
	if(label in selected_labels):
		color, marker = selected_labels[label]
		plt.scatter(expr, lfc, color=color, marker=marker, s = 72, label = label)
	else:
		plt.scatter(expr, lfc, s = 2)
		

legend = plt.legend(scatterpoints = 1, loc=(0.7, 0.7), prop={'size':10}, frameon = False)
#legend.get_frame().set_edgecolor('white')
plt.ylabel('log2(KO/WT)')
plt.xlabel('log2(RPKM in WT)')
if(args.output):
	plt.savefig(args.output);
else:
	plt.show()
	