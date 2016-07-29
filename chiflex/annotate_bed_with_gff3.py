#! /usr/bin/python
'''assignes genomic features to intervals they overlap with''' 
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool, Interval

from nrlbio.pybedtools_extension import construct_gff_interval, gff2bed, bed2gff
from nrlbio.genome_system import gff3_to_genes


parser = argparse.ArgumentParser(description='assignes genomic features to intervals they overlap with');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file to be annotated");
parser.add_argument('--gff3', nargs = '?', required = True, type = str, help = "path to annotation, gff3 file");
#parser.add_argument('--genes', nargs = '?', required = True, type = str, help = "path to annotation, yml file. Ignored if [--gff3] is set");

parser.add_argument('-f', '--fraction', nargs = '?', default = 0.5, type = float, help = "if fraction of bed interval covered by an annotation feature is more than set value, then the feature is reported for the interval. Corresponds to \'-f\' option of intersectBed");
parser.add_argument('--feature', nargs = '?', default = 'interval', type = str, help = "feature name for all entries in the annotated file");
args = parser.parse_args();

genes = dict([(x.name, x) for x in gff3_to_genes(BedTool(args.gff3))])
gene_intervals = [gff2bed(x.gene) for x in genes.values()]


bed = BedTool(args.path);
offset  = bed.field_count();
intersection = bed.intersect(b=gene_intervals, s=True, wao=True, f=args.fraction);
annotations = defaultdict(list)

for interval in intersection:
	name = interval[offset+3]
	if(name == '.'):
		annotations[interval.name].append(('intergenic', '', '', '', ''))
	else:
		gene = genes.get(name.split(":")[1]);
		regulation, transcription, biotypes  = gene.annotate_interval(interval);
		annotations[interval.name].append((','.join(biotypes), ','.join(transcription), ','.join(regulation), gene.name, gene.gene_symbol))
	

if(bed.file_type=='bed'):
	for interval in bed:
		annotation = annotations[interval.name]
		
		regulation = ":".join([x[0] for x in annotation])
		transcription = ":".join([x[1] for x in annotation])
		biotypes = ":".join([x[2] for x in annotation])
		names = ":".join([x[3] for x in annotation])
		gene_symbols = ":".join([x[4] for x in annotation])
		
		attrs =  (('regulation', regulation), ('transcription', transcription), ('biotypes', biotypes), ('names', names), ('gene_symbols', gene_symbols))
		sys.stdout.write(str(bed2gff(interval, feature='an', attrs=attrs)))
		
elif(bed.file_type =='gff'):
	for interval in bed:
		annotation = annotations[interval.name]
		
		regulation = ":".join([x[0] for x in annotation])
		transcription = ":".join([x[1] for x in annotation])
		biotypes = ":".join([x[2] for x in annotation])
		names = ":".join([x[3] for x in annotation])
		gene_symbols = ":".join([x[4] for x in annotation])
		
		for k, v in (('regulation', regulation), ('transcription', transcription), ('biotypes', biotypes), ('names', names), ('gene_symbols', gene_symbols)):
			interval.attrs[k]=v;
		sys.stdout.write(str(interval))
else:
	sys.exit('Only bed or gff files can be annotated yet:(\n');




