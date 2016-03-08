#! /usr/bin/python
'''Selects for targets with a destructive potential''' 
import argparse
import sys;
import math;
import copy

from Bio import SeqIO;
import jinja2

from nrlbio.generators import generator_doublebed
from nrlbio.rnahybrid import get_rnahybrid, get_shuffled_rnahybrid
from nrlbio.html import get_link, get_tarbase_link
from nrlbio.pybedtools_extension import construct_gff_interval

system_choices = ['hg19', 'hg38', 'mm9', 'mm10', 'ce6', 'circ']
sys2rhsys = {'hg19': '3utr_human', 'hg38': '3utr_human', 'mm9': '3utr_human', 'mm10': '3utr_human', 'ce6': '3utr_worm', 'ce10': '3utr_worm', 'circ': '3utr_human'}

parser = argparse.ArgumentParser(description='Selects for targets with a destructive potential');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mirna:target interactions to select from, double gff/bed. Sequences have to be assigned to each interaction");
parser.add_argument('--system', nargs = '?', required = True, choices = system_choices, type = str, help = "Genome system. Can be set to %s" % "|".join(system_choices))
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the miRNAs, fasta format");
parser.add_argument('-o', '--output', nargs = '?', required = True, type = str, help = "Path to the output report");
parser.add_argument('--internal', nargs = '?', const = True, default = False, type = bool, help = "If set, links will point to internal UCSC browser")
parser.add_argument('--template', nargs = '?', default = "interaction_destructive.html", type = str, help = "path to jinja2/django template")
parser.add_argument('--title', nargs = '?', default = 'Destructive Interactions', type = str, help = "title of html report");
parser.add_argument('--pval_cutoff', nargs = '?', default = 0.1, type = float, help = "Max p value required for additional(split) interactions to be output");
parser.add_argument('--best_only', nargs = '?', default = 0, type = int, help = "If set, only the top [best_only] interactions will be output");
args = parser.parse_args();

rhsys = sys2rhsys[args.system]

mirdict = {};
for seqrecord in SeqIO.parse(args.mir, "fasta"):
	mirdict[seqrecord.id] = str(seqrecord.seq.upper())
	
	
class Interaction(object):
	def __init__(self, interval, basepairing, system, internal):
		self.interval = interval;
		self.basepairing = basepairing;
		self.targetpairing = ("   target  %s" % basepairing[0], "   5\'      %s" % basepairing[1])
		self.mirnapairing = ("           %s  miRNA" % basepairing[2], "           %s  5\'" % basepairing[3])
		self.ilink = get_link(interval, system, internal = internal)
		self.mlink = get_tarbase_link(interval.name)	


def basepairing2seq(basepairing):
	s = []
	for u, a in zip(basepairing[0], basepairing[1]):
		if(a!=' '):
			s.append(a);
		if(u!=' '):
			s.append(u);
	return ''.join(s).replace('U', 'T')


def split_interval(interval, seq, length_limit=20):
	start = interval.attrs['seq'].find(seq)
	end = start+len(seq);
	add = (length_limit - len(seq))/2
	if(add>0):
		start = max(0, start-add);
		end = min(end+add, len(interval));
		seq = interval.attrs['seq'][start:end]
	
	if(interval.strand == '-'):
		lboundaries = (interval.end-start, interval.end, interval.attrs['seq'][:start])
		rboundaries = (interval.start, interval.end-end, interval.attrs['seq'][end:])
		interval.start = interval.end - end
		interval.end = interval.end - start
	else:
		lboundaries = (interval.start, start, interval.attrs['seq'][:start])
		rboundaries = (interval.start+end, interval.end, interval.attrs['seq'][end:])
		interval.end = interval.start + end;
		interval.start = interval.start + start 

		
	interval.attrs['seq'] = seq;
	
	r = [];
	for b in filter(lambda x: x[1]-x[0]>=length_limit, (lboundaries, rboundaries)):
		r.append(construct_gff_interval(interval.chrom, b[0], b[1], 'ai', score='0', strand=interval.strand, attrs=[ ('ID', interval.name), ('seq', b[2]), ('n_uniq', (int(interval.attrs['n_uniq'])+1)/2), ('split', 'yes') ] ))
		
	return r;

		
		
def pattern2score(pattern):
	return sum(pattern[1:9]) - sum(pattern[9:11]) + sum(pattern[11:19])

def interval2score(interval, energy, pattern):
	return pattern2score(pattern) - energy + math.log(int(interval.attrs['n_uniq']));


def get_interaction(interval, mirseq, name, gsystem, pval_cutoff):
	interval.name = name
	energy, pattern, basepairing, pval = get_rnahybrid(interval.attrs['seq'], mirseq, system = rhsys, extended=True);
	if(pval>pval_cutoff):
		return None, []
	else:
		interval.attrs['dscore'] = "%1.2f" % interval2score(interval, energy, pattern)
		interval.attrs['pval'] = "%1.5f" % pval
		interval.attrs['energy'] = str(energy);
		interval.attrs['pattern'] = ",".join([str(x) for x in pattern]);
		interval.attrs['mseq'] = mirseq;
		interval.attrs['split'] = interval.attrs.get('split', 'no')
		seq = basepairing2seq(basepairing)
		aintervals = split_interval(interval, seq, length_limit=20);
		sys.stdout.write(str(interval))
		return Interaction(interval, basepairing, gsystem, args.internal), aintervals


interactions = [];
for count, (i1, i2) in enumerate(generator_doublebed(args.path)):
	name = i1.chrom
	mirseq = mirdict[i1.chrom]
	#ti = copy.copy(i2)
	interaction, aintervals = get_interaction(i2, mirseq, name, args.system, 1)
	interactions.append(interaction)
	while(aintervals):
		newaints = [];
		for ainterval in aintervals:
			interaction, aints = get_interaction(ainterval, mirseq, name, args.system, args.pval_cutoff)
			newaints.extend(aints)
			if(interaction):
				interactions.append(interaction)
		aintervals = newaints;
		
	if(count and count % 10000 == 0):
		sys.stderr.write("%d interactions are processed\n" % count);
			


	

interactions.sort(key=lambda x: float(x.interval.attrs['dscore']), reverse=True)
if(args.best_only):
	interactions = interactions[:args.best_only]


	
environment = jinja2.Environment(loader=jinja2.PackageLoader('nrlbio', 'templates'))
template = environment.get_template(args.template)
with open(args.output, 'w') as ofile:
	ofile.write(template.render({"interactions": interactions, 'title': args.title}));