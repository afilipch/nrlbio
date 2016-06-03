#! /usr/bin/python
'''Annotates circula RNA with an ENSEMBL annotation'''


import sys;
import argparse
from collections import defaultdict

from pybedtools import BedTool, Interval

#from nrlbio.pybedtools_extension import construct_gff_interval(chrom, start, stop, feature, score='0', strand='.', source='un', frame='.', attrs=[])


parser = argparse.ArgumentParser(description='Annotates circula RNA with an ENSEMBL annotation');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the circles, bed/gff format");
parser.add_argument('--gff3', nargs = '?', required = True, type = str, help = "Path to the genome system annotation, gff3 format");
parser.add_argument('--convert_name', nargs = '?', default = False, const = True, type = bool, help = "If set, ENSEMBL chromosme names will be converted to the UCSC format");
args = parser.parse_args();

def get_exonic_structure(exons, tnames):
	res = [];
	temp = [];
	for tname in tnames:
		for exon in exons:
			if(tname in exon[4]):
				temp.append("-".join([exon[1], exon[2]]))
		res.append(tuple(temp));
		temp = [];
		
	return list(set(res));




def get_exons(interval, exons):
	'''Outputs all exons which may form a structure of the given interval(circle)
	
	interval pybedtools.Interval: genomic interval representing a circle (boundaries of the transcript)
	exons list: each element is a tuple representing genomic interval of an exon
	'''
	
	exons = sorted(list(set(exons)), key = lambda x: x[1])
	exons = filter(lambda x: int(x[1]) >= interval.start and int(x[2]) <= interval.end, exons)
	
	exonic_boundaries = False
	same_gene = False
	same_transcript = False
	texons = []
	tnames = []
	
	
	if(exons):
		nexons = []
		for exon in exons:
			nexon = list(exon)
			nexon[3] = exon[3].split(",")
			nexon[4] = exon[4].split(",")
			nexons.append(nexon)
			
		ef, el = nexons[0], nexons[-1]
		if(interval.start == int(ef[1]) and  interval.end == int(el[2])):
			exonic_boundaries = True
			gnames = set(ef[3]) & set(el[3])
			if(gnames):
				same_gene = True	
			tnames = set(ef[4]) & set(el[4])
			if(tnames):
				same_transcript = True
				texons = get_exonic_structure(nexons, tnames)

		
		
		
	#for exon in exons:
		#print "\t".join(exon);
		
	#print 
	#print "%s\t%d\t%d\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.end, interval.name, '0', interval.strand) 
	#print exonic_boundaries, same_gene, same_transcript, tnames, "|".join([",".join(x) for x in texons])
	#print
	#print "_"*130
	
	interval.attrs['exonic_boundaries'] = str(exonic_boundaries)
	interval.attrs['same_gene'] = str(exonic_boundaries)
	interval.attrs['same_transcript'] = str(exonic_boundaries)
	interval.attrs['transcripts'] = "|".join(tnames)
	interval.attrs['exonic_structure'] = "|".join([",".join(x) for x in texons])
	print "\t".join(str(interval).split("\t")[:9])
		
	return None




#get exons from annotation gff3 file
exon2names = defaultdict(list);

for interval in BedTool(args.gff3):
	if('ID' in interval.attrs and interval.attrs['ID'].split(':')[0] == 'gene'):
		gname = interval.attrs['gene_id']
	elif('ID' in interval.attrs and interval.attrs['ID'].split(':')[0] == 'transcript'):
		tname = interval.attrs['transcript_id']
		
	elif(interval[2] == 'exon'):
		if(args.convert_name):
			if(interval.chrom == 'MT'):
				interval.chrom = 'chrM'
			else:
				interval.chrom = "chr%s" % interval.chrom
		exon2names[(interval.chrom, interval.start, interval.end, interval.strand)].append((gname, tname))
		

exons = []
for gi, names in exon2names.iteritems():
	exons.append(Interval(gi[0], gi[1], gi[2], name=','.join(set([x[0] for x in names])), score=','.join(set([x[1] for x in names])), strand=gi[3]))

	
#Get an intersection between circles and exons
bed = BedTool(args.path);
bed.sort()
offset  = bed.field_count();
intersection = bed.intersect(b=exons, s=True, wao=True);

curname = ''
cexons = []
for interval in intersection:
	if(curname == interval.name):
		cexons.append(tuple(interval[offset:offset+6]))
	else:
		if(curname):
			get_exons(cinterval, cexons)
		cinterval = interval
		cexons = [tuple(interval[offset:offset+6])]
		curname = interval.name
else:
	get_exons(cinterval, cexons)



