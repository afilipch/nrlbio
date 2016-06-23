#! /usr/bin/python
'''Selects for targets with a destructive potential''' 
import argparse
import sys;
import math;
import copy

from Bio import SeqIO;
import jinja2
from pybedtools import BedTool

from nrlbio.mirna import slicing_score
from nrlbio.rnahybrid import get_rnahybrid
from nrlbio.html import get_link, get_tarbase_link


system_choices = ['hg19', 'hg38', 'mm9', 'mm10', 'ce6', 'ce11', 'circ']
sys2rhsys = {'hg19': '3utr_human', 'hg38': '3utr_human', 'mm9': '3utr_human', 'mm10': '3utr_human', 'ce6': '3utr_worm', 'ce10': '3utr_worm', 'ce11': '3utr_worm', 'circ': '3utr_human'}

parser = argparse.ArgumentParser(description='Selects for targets with a destructive potential');

parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mirna:target interactions to select from, double gff/bed. Sequences have to be assigned to each interaction");
parser.add_argument('--system', nargs = '?', required = True, choices = system_choices, type = str, help = "Genome system. Can be set to %s" % "|".join(system_choices))
parser.add_argument('-o', '--output', nargs = '?', required = True, type = str, help = "Path to the output report");

parser.add_argument('--internal', nargs = '?', const = True, default = False, type = bool, help = "If set, links will point to internal UCSC browser")
parser.add_argument('--template', nargs = '?', default = "interaction_destructive.html", type = str, help = "path to jinja2/django template")
parser.add_argument('--title', nargs = '?', default = 'Destructive Interactions', type = str, help = "title of html report");
parser.add_argument('--best_only', nargs = '?', default = 0, type = int, help = "If set, only the top [best_only] interactions will be output");
parser.add_argument('--conservation', nargs = '?', default = False, const = True, type = bool, help = "If set, conservation score will be taken into account");
args = parser.parse_args();


rhsys = sys2rhsys[args.system]

class Interaction(object):
	def __init__(self, interval, basepairing, system, internal):
		self.interval = interval;
		self.basepairing = basepairing;
		self.targetpairing = ("   target  %s" % basepairing[0], "   5\'      %s" % basepairing[1])
		self.mirnapairing = ("           %s  miRNA" % basepairing[2], "           %s  5\'" % basepairing[3])
		self.ilink = get_link(interval, system, internal = internal)
		self.mlink = get_tarbase_link(interval.attrs['mirid'])
		
	def __cmp__(self, other):
		for key in ('destructive_score', 'cons_dscore'):
			a = cmp(float(self.interval.attrs.get(key, 0)), float(other.interval.attrs.get(key, 0)))
			if(a):
				return a
		return a	



if(args.conservation):
	def interval2score(interval):
		cons_dscores, cons_bscores, cons_bls = [interval.attrs[x].split(',') for x in ('cons_dscores', 'cons_bscores', 'cons_bls')]
		cons_score = 0;
		for dscore, bscore, bl in zip(cons_dscores, cons_bscores, cons_bls):
			cons_score += float(dscore)*float(bscore)*math.log(float(bl) + 4, 4)
		return float(interval.attrs['dscore']) + cons_score/math.log(2+len(dscore), 4)
else:
	def interval2score(interval):
		return float(interval.attrs['dscore'])


interactions = [];
for count, interval in enumerate(BedTool(args.path)):
	
	energy, pattern, basepairing, pval, pos = get_rnahybrid(interval.attrs['seq'], interval.attrs['mseq'], system = rhsys, extended=True);
	
	interval.attrs['destructive_score'] = "%1.2f" % interval2score(interval) 
	interval.attrs['slicing_score'] = "%1.2f" % slicing_score(pattern)
	interval.attrs['energy'] = "%1.2f" % energy

	
	interactions.append(Interaction(interval, basepairing, args.system, args.internal))
	sys.stdout.write(str(interval))
	
	if(count and count % 10000 == 0):
		sys.stderr.write("%d interactions are processed\n" % count);
		




interactions.sort(reverse=True)
#interactions.sort(key=lambda x: float(x.interval.attrs['dscore']), reverse=True)
if(args.best_only):
	interactions = interactions[:args.best_only]


environment = jinja2.Environment(loader=jinja2.PackageLoader('nrlbio', 'templates'))
template = environment.get_template(args.template)
with open(args.output, 'w') as ofile:
	ofile.write(template.render({"interactions": interactions, 'title': args.title}));
	
	
	
	
	
	