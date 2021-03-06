# /usr/bin/python
'''Draws a piechart of biotypes distribuiton of provided genomic intervals'''
import sys
import argparse
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.pyplot_extension import pie


parser = argparse.ArgumentParser(description='Draws a piechart of biotypes distribuiton of provided genomic intervals');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated genomic intervals, gff file");
parser.add_argument('-hi', '--hierarchical', nargs = '?', default=False, const = True, type = bool, help = "If set, hierarchical assignemt of confilicting annotations will be applied, uniform assignment otherwise");
parser.add_argument('-r', '--reads', nargs = '?', default=False, const = True, type = bool, help = "If set, number of chimeras, not number of interactions will be counted");
parser.add_argument('--output', nargs = '?', default='biotypes.png', type = str, help = "Path to the output");
parser.add_argument('--title', nargs = '?', type = str, help = "Title for a plot");
args = parser.parse_args();


def collapse_hierarchical(btypes, ttypes, score):
	types  = set()
	for blist,tlist in zip(btypes, ttypes):
		if('exon' in tlist):
			types.update(blist);
	
	for region in ['utr3', 'utr5', 'cds']:
		if(region in types):
			return {region:score}
			
	return {'intron': 0}
			


def collapse_uniform(btypes, ttypes, score):
	blists  = []
	for blist,tlist in zip(btypes, ttypes):
		if('exon' in tlist):
			blists.append(blist);
	
	ans = defaultdict(float)
	gnorm = float(len(blists));
	for blist in blists:
		lnorm = gnorm*len(blist);
		for b in blist:
			ans[b] += score/lnorm;
			
	return ans;
			


if(args.hierarchical):
	collapse = collapse_hierarchical
else:
	collapse = collapse_uniform
			

biotypes = defaultdict(float);
for interval in BedTool(args.path):
	coding = [] 
	for p, c in enumerate(interval.attrs['regulation'].split(':')):
		if('protein_coding' in c.split(',')):
			coding.append(p);

			
	btypes = interval.attrs['biotypes'].split(':')
	btypes = [btypes[x].split(',') for x in coding if btypes[x]]
	ttypes = interval.attrs['transcription'].split(':')
	ttypes = [ttypes[x].split(',') for x in coding if ttypes[x]]


	if(args.reads):
		score = float(interval.attrs['n_uniq'])
	else:
		score = 1.0
		
	if(btypes):
		for k, v in collapse(btypes, ttypes, score).items():
			if(v):
				biotypes[k] += v;
		

#print sum(biotypes.values())
		
for k,v in biotypes.items():
	print "%s\t%d" % (k,v)

	
pie(biotypes, top=10, min_fraction=0.03, title=args.title, labelsize='medium', labelva='top', output=args.output, explode=None, labels=None, colors=('yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'azure', 'seashell', 'darkorchid', 'chartreuse'), pctdistance=0.8, shadow=False, labeldistance=1.15, startangle=0, radius=None, counterclock=True, wedgeprops=None, textprops={'fontsize': 'medium'}, frame=(0.12, 0.12, 0.76, 0.76));