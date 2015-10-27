#! /usr/bin/python
'''assignes genomic features to intervals they overlap with''' 
import argparse
import sys;
from collections import defaultdict, namedtuple
from itertools import groupby 
from operator import itemgetter, attrgetter

from pybedtools import BedTool, Interval

from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='assignes genomic features to intervals they overlap with');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file to be annotated");
parser.add_argument('-a', '--annotation', nargs = '?', required = True, type = str, help = "path to annotation gff file");
parser.add_argument('-f', '--fraction', nargs = '?', default = 0.5, type = float, help = "if fraction of bed interval covered by an annotation feature is more than set value, then the feature is reported for the interval. Corresponds to \'-f\' option of intersectBed");
parser.add_argument('--feature', nargs = '?', default = 'interval', type = str, help = "feature name for all entries in the annotated file");
args = parser.parse_args();

def get_name(ins, offset):
	d = dict([x.split('=') for x in ins[offset+8].split(';')]);
	return d['Name'];


genes = [];
for interval in BedTool(args.annotation):
	gid = interval.attrs.get('ID');
	if(gid):
		ftype, name = gid.split(":");
		if(ftype == 'gene'):
			interval.chrom = "chr%s" % interval.chrom
			genes.append(interval)
			#sys.stdout.write(str(interval));
			#print gid, interval.name
			#print
			
			
genes = BedTool(genes);
bed = BedTool(args.path);
offset  = bed.field_count();

intersection = bed.intersect(b=genes, s=True, wo=True, f=args.fraction);
curname = intersection[0].name;
gene_names = [get_name(intersection[0], offset)];

for ins in intersection[1:]:
	if(ins.name == curname):
		gene_names.append(get_name(ins, offset));
	else:
		ins.attrs['genes'] = ",".join(gene_names)
		sys.stdout.write(str(ins));
		curname = ins.name;
		gene_names = [get_name(ins, offset)];
else:
	ins.attrs['genes'] = ",".join(gene_names)
	sys.stdout.write(str(ins));


