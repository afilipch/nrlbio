# /usr/bin/python
'''Draws a piechart of regulation distribuiton of provided genomic intervals'''
import sys
import argparse
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.pyplot_extension import pie


parser = argparse.ArgumentParser(description='Draws a piechart of regulation distribuiton of provided genomic intervals');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated genomic intervals, gff file");
parser.add_argument('-hi', '--hierarchical', nargs = '?', default=False, const = True, type = bool, help = "If set, hierarchical assignemt of confilicting annotations will be applied, uniform assignment otherwise");
parser.add_argument('-r', '--reads', nargs = '?', default=False, const = True, type = bool, help = "If set, number of chimeras, not number of interactions will be counted");
parser.add_argument('--output', nargs = '?', default='regulation.png', type = str, help = "Path to the output");
parser.add_argument('--title', nargs = '?', type = str, help = "Title for a plot");
args = parser.parse_args();


ORDER = ['protein_coding', 'miRNA', 'scaRNA', 'snoRNA', 'rRNA', 'Mt_tRNA', 'Mt_rRNA', 'snRNA', 'lincRNA', 'intron', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene',  'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene', 'polymorphic_pseudogene', 'processed_transcript', '3prime_overlapping_ncrna', 'misc_RNA', 'non_coding', 'nonsense_mediated_decay', 'non_stop_decay', 'retained_intron', 'sense_overlapping',  'antisense',  'sense_intronic', 'TEC',  'intergenic'] 

#NONEXONIC = ['nonsense_mediated_decay', 'non_stop_decay', 'retained_intron', 'sense_overlapping',  'antisense',  'sense_intronic', 'TEC',  'intergenic'];

def collapse_hierarchical(regulation, transcription, score):
	variants = set()
	
	for reg, tr in zip(regulation, transcription):
		if('exon' in tr):
			variants.update(reg)
		elif('intron' in tr):
			variants.add('intron')
	
	if(not variants):
		for reg in regulation:
			variants.update(r);
					
	if(not variants):
		return {'unknown':score}
	
	
	for o in ORDER:
		if(o in variants):
			return {o:score}




def collapse_uniform(regulation, score):
	ans = defaultdict(float)
	norm = len(regulation)
	
	for reg in regulation:
		tnorm = norm*len(reg)
		for r in reg:
			ans[r] += (1.0*score)/tnorm
			
	return ans


if(args.hierarchical):
	collapse = collapse_hierarchical
else:
	collapse = collapse_uniform
	

types = defaultdict(float);
for interval in BedTool(args.path):
	regulation = [x.split(',') for x in interval.attrs['regulation'].split(':') if x];
	#variants = set()
	#for reg in regulation:
		#variants.update(reg)
	
	transcription= [x.split(',') for x in interval.attrs['transcription'].split(':') if x];
	#ttypes = set();
	#for t in transcription:
		#ttypes.update(t);
	
	#if(len(ttypes)==1 and list(ttypes)[0] == 'intron'):
		#regulation.append(['intron'])
	
	
	#print regulation 
	#print 
	#print biotypes
	#print 
	#print "-"*140
	
	if(args.reads):
		score = float(interval.attrs['n_uniq'])
	else:
		score = 1.0
		
	for k,v in collapse(regulation, transcription, score).items():
		types[k] += v;
		
			
for k in ORDER:
	print "%s\t%d" % (k, types[k])
	
	
pie(types, top=10, min_fraction=0.03, title=args.title, labelsize='medium', labelva='top', output=args.output, explode=None, labels=None, colors=('yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'azure', 'seashell', 'darkorchid', 'chartreuse'), pctdistance=0.8, shadow=False, labeldistance=1.15, startangle=0, radius=None, counterclock=True, wedgeprops=None, textprops={'fontsize': 'medium'}, frame=(0.12, 0.12, 0.76, 0.76));	
			