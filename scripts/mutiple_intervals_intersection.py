#! /usr/bin/python
'''Detects how many intervals are shared between all possible pairs between provided bed/gff files''' 
import sys;
import os;
import argparse;

from pybedtools import BedTool;
from itertools import combinations


parser = argparse.ArgumentParser(description='Detects how many intervals are shared between all possible pairs between provided bed/gff files');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "paths to the bed/gff files to intersect");
parser.add_argument('--names', nargs = '+', type = str, help = "Names corresponding to the bed/gff files provided. If set, the number of names has to be equal to the number of bed/gff files. If not set, bed/gff files base names will be used instead");
parser.add_argument('-f', '--fraction', nargs = '?', default = 0.95, type = float, help = "fraction(in nucleotides) that intervals have to reciprocally overlap one on another to be accounted shared between two bed/gff files")
parser.add_argument('--drawdir', nargs = '?', type = str, help = "Path to the folder to write venn diagrams in. If not set, plots will not be produced");
args = parser.parse_args()

if(args.names):
	if(len(args.names) == len(args.path)):
		names = args.names;
	else:
		sys.stderr.write("WARNING: The number of names provided has to be equal to the number of bed/gff files. Bed/gff files base names will be used instead\n")
		names = [os.path.basename(x) for x in args.path]
else:
	names = [os.path.basename(x) for x in args.path]
	
	
results = [];	
bedtools = [BedTool(x) for x in args.path];
for i1, i2 in combinations(list(range(len(bedtools))), 2):
	intersection = bedtools[i1].intersect(bedtools[i2], u=True, s=True, f=args.fraction, F=args.fraction);
	results.append((names[i1], names[i2], len(bedtools[i1]), len(bedtools[i2]), len(intersection)));
	
for r in results:
	print "\t".join([str(x) for x in r])






if(args.drawdir):
	from matplotlib import pyplot as plt
	from matplotlib_venn import venn2, venn2_circles
	
	colors = ('yellow', 'blue');

	try:
		os.mkdir(args.drawdir)
	except:
		sys.stderr.write("directory %s already exists, plots will be (over)written there\n" % os.path.abspath(args.drawdir));
		
	def get_plot(result, colors, outdir):
		subsets = (result[2]-result[4], result[3]-result[4], result[4]);
		labels = result[:2];
		title = "_VS_".join(labels);
		v = venn2(subsets=subsets, set_labels=labels, set_colors=colors)
		#v.get_patch_by_id('10').set_color(colors[0])
		#v.get_patch_by_id('01').set_color(colors[1])
		plt.savefig(os.path.join(outdir, title))
		plt.close()
		
		
	for result in results:
		get_plot(result, colors, args.drawdir)
		
	
		
		
		
		
		





## Subset labels
#v.get_label_by_id('10').set_text('A but not B')
#v.get_label_by_id('01').set_text('B but not A')
#v.get_label_by_id('11').set_text('A and B')

## Subset colors
#v.get_patch_by_id('10').set_color('c')
#v.get_patch_by_id('01').set_color('#993333')
#v.get_patch_by_id('11').set_color('blue')

## Subset alphas
#v.get_patch_by_id('10').set_alpha(0.4)
#v.get_patch_by_id('01').set_alpha(1.0)
#v.get_patch_by_id('11').set_alpha(0.7)

## Border styles
#c = venn2_circles(subsets=s, linestyle='solid')
#c[0].set_ls('dashed')  # Line style
#c[0].set_lw(2.0)       # Line width

		
		
		
		
		
		
		
		
		
		
		
		
		
		