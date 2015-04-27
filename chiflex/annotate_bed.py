#! /usr/lib/python
'''assignes genomic features to intervals they overlap with''' 
import argparse
import sys;
from collections import defaultdict, namedtuple
from itertools import groupby 
from operator import itemgetter, attrgetter

from pybedtools import BedTool, Interval, create_interval_from_list

from nrlbio.pybedtools_extension import construct_gff_interval

#a = [[1, 7], [4, 6], [1, 3], [4, 0], [0, 0]]
#for i in groupby(sorted(a, key=itemgetter(0)), key=itemgetter(0)):
	#for x in i[1]:
		#print i[0], x;
#sys.exit()

parser = argparse.ArgumentParser(description='assignes genomic features to intervals they overlap with');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file to be annotated");
parser.add_argument('-a', '--annotation', nargs = '?', required = True, type = str, help = "path to annotation gff file");
parser.add_argument('-f', '--fraction', nargs = '?', default = 0.5, type = int, help = "if fraction of bed interval covered by an annotation feature is more than set value, then the feature is reported for the interval. Corresponds to \'-f\' option of intersectBed");
parser.add_argument('--feature', nargs = '?', default = 'interval', type = str, help = "feature name for all entries in the annotated file");
args = parser.parse_args();


Annotation = namedtuple('Annotation', 'type, transcription, regulation, gene_id, score')

def intersection2annotation(interval, offset):
	if(interval[offset] != "."):
		g = create_interval_from_list(interval[offset:]);
		return Annotation(g[2], g.attrs['transcription'], g.attrs['regulation'], g.attrs['gene_id'], float(g[-1]))
	else:
		return Annotation("intergenic", "", "", "", 0)
		
		
def collapse_annotation(annotations, interval, best_only=True):
	if(best_only):
		bs = max([x.score for x in annotations]);
		filtered = filter(lambda x: x.score==bs, annotations)
	else:
		filtered = annotations
		
	types = [];
	regulation = [];
	transcription = [];
	processed = list(sorted(list(set(filtered)), key = attrgetter('type', 'regulation', 'transcription')))
	
	for c1, (k1, l1) in enumerate(groupby(processed, key = attrgetter('type'))):
		types.append(k1);
		regulation.append([])
		transcription.append([])
		for c2, (k2, l2) in enumerate(groupby(l1, key = attrgetter('regulation'))):
			regulation[c1].append(k2)
			transcription[c1].append([]);
			for a in l2:
				transcription[c1][c2].append(a.transcription);
				
	ts = ",".join(types)
	rn = ",".join(["|".join(x) for x in regulation])
	tn = ",".join(["|".join(["%".join(y) for y in x]) for x in transcription])
	gene_ids = ",".join(set([x.gene_id for x in filtered]))
	
	ad = interval.attrs;
	for k, v in (('type', ts), ('regulation', rn), ('transcription', tn), ('gene_name', gene_ids), ('ID', interval.name)):
		ad[k] = v;
	
	return construct_gff_interval(interval.chrom, interval.start, interval.stop, feature=args.feature, score=interval.score, strand=interval.strand, source='annotate_bed.py', frame='.',
	attrs= ad.items())
	
	
bed = BedTool(args.path).sort();
offset  = bed.field_count();
gff_annotation = BedTool(args.annotation).sort();

intersection = bed.intersect(gff_annotation, s=True, wao=True, f=args.fraction, sorted = True, nonamecheck=True)
previous = intersection[0];
curname = previous.name;
local_annotation = [intersection2annotation(previous, offset)];

for i in intersection:
	if(i.name == curname):
		local_annotation.append(intersection2annotation(i, offset));
	else:
		sys.stdout.write(str(collapse_annotation(local_annotation, previous)));
		curname = i.name;
		previous = i;
		local_annotation = [intersection2annotation(i, offset)];
else:
	sys.stdout.write(str(collapse_annotation(local_annotation, previous)))