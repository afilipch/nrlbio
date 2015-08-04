#! /usr/lib/python
'''evaluates chiflex performance on artificially generated reads'''
import argparse
import os;
import sys;
from collections import defaultdict, namedtuple, Counter
from itertools import chain

import pysam;
import yaml
from Bio import SeqIO;

from nrlbio.numerictools import overlap as get_overlap
from nrlbio.generators import generator_doublesam
from nrlbio.samlib import sort_segments_qstart


parser = argparse.ArgumentParser(description='evaluates chiflex performance on artificially generated reads');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to initial reads. Fasta file");
parser.add_argument('-su', '--single_unique', nargs = '?', required = True, type = str, help = "Path to the unique single hits. Sam/bam format");
parser.add_argument('-sc', '--single_control', nargs = '?', required = True, type = str, help = "Path to the control single hits. Sam/bam format");
parser.add_argument('-cu', '--chimera_unique', nargs = '?', required = True, type = str, help = "Path to the unique chimera hits. Sam/bam format");
parser.add_argument('-cc', '--chimera_control', nargs = '?', required = True, type = str, help = "Path to the control chimera hits. Sam/bam format");
parser.add_argument('-sn', '--single_nonunique', nargs = '?', required = True, type = str, help = "Path to the nonunique single hits. Sam/bam format");
parser.add_argument('-cn', '--chimera_nonunique', nargs = '?', required = True, type = str, help = "Path to the nonunique chimera hits. Sam/bam format");
parser.add_argument('--coverage_cutoff', nargs = '?', default = 0.8, type = float, help = "if fraction of read covered by mapping is more than [coverage_cutoff], then such a mapping is reported to be exact");
parser.add_argument('-o', '--outdir', nargs = '?', default = 'evaluation_data', type = str, help = "Name of the output directory. If directory exists, files will be (re)written there")
args = parser.parse_args();


Chimera = namedtuple('Chimera', 'typename chrom1 strand1 start1 stop1 ttype1 chrom2 strand2 start2 stop2 ttype2 left right gap conv1 conv2')
Single = namedtuple('Single', 'typename chrom strand start stop ttype left right conv')
Shuffled = namedtuple('Shuffled', 'typename chrom strand start stop ttype length')
Random = namedtuple('Random', 'typename length ttype');

#get initial reads
#_____________________________________________________________________________________________________________________________
def name2chimera(region, ttype, gaps, convnum):
	chrom1, strand1, start1, stop1, chrom2, strand2, start2, stop2= chain(*[x.split(":") for x in region.split("&")]);
	start1, stop1, start2, stop2 = [int(x) for x in (start1, stop1, start2, stop2)]
	conv1, conv2 = [int(x) for x in convnum.split(":")]
	
	ttype1, ttype2 = ttype.split(":");
	left, right, gap = [int(x) for x in gaps.split(":")]
	
	return Chimera('chimera', chrom1, strand1, start1, stop1, ttype1, chrom2, strand2, start2, stop2, ttype2, left, right, gap, conv1, conv2)
	
	
def name2single(region, ttype, gaps, convnum):
	chrom, strand, start, stop = region.split(":");
	start, stop = int(start), int(stop)
	conv = int(convnum)

	left, right = [int(x) for x in gaps.split(":")]
	
	return Single('single', chrom, strand, start, stop, ttype, left, right, conv)
	
	
def name2shuffled(region, ttype, gaps):
	chrom, strand, start, stop = region.split(":");
	start, stop = int(start), int(stop)

	length = int(gaps)
	
	return Shuffled('shuffled', chrom, strand, start, stop, ttype, length)	
	
	
def name2random(gaps):
	length = int(gaps)	
	return Random('random', length, 'random')		
	

def get_source(name):
	number, rtype, region, ttype, gaps, convnum = name.split("|")
	if(rtype == "chimera"):
		return name2chimera(region, ttype, gaps, convnum)
	elif(rtype == "single"):
		return name2single(region, ttype, gaps, convnum)
	elif(rtype == 'shuffled'):
		return name2shuffled(region, ttype, gaps)
	elif(rtype == 'random'):
		return name2random(gaps)
	else:
		sys.exit('Unrecognized read type. Has to be chimera|single|shuffled|random')
initial = defaultdict(int);



reads = [];
for i, seqrecord in enumerate(SeqIO.parse(args.path, 'fasta')):
	reads.append(get_source(seqrecord.name));
	
	
	
	
	
	
#_____________________________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________________________
#get mapped reads
#_____________________________________________________________________________________________________________________________	
#_____________________________________________________________________________________________________________________________	
ComparisonSingle = namedtuple('ComparisonSingle', 'num length read_type mapped_type overlap');
ComparisonChimera = namedtuple('ComparisonChimera', 'num length1 length2 read_type mapped_type overlap1 overlap2');
ComparisonControl = namedtuple('ComparisonControl', 'num length read_type mapped_type AS');

strand_conv = {'-': True, '+': False}

comparison_dict = defaultdict(list);	





#_____________________________________________________________________________________________________________________________
#unique single reads	
#_____________________________________________________________________________________________________________________________
def single_unique_comparison(rnum, read, segment, segchrom):
	mapped_type = 'single'
	if(read.typename == 'single'):
		if(read.chrom == segchrom and strand_conv[read.strand] == segment.is_reverse):
			s, e = get_overlap((read.start, read.stop), (segment.reference_start, segment.reference_end));
			overlap = max(e-s, 0) 
			return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type, overlap);
		else:
			return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type, 0);
			
	elif(read.typename == 'chimera'):
		if(read.chrom1 == segchrom and strand_conv[read.strand1] == segment.is_reverse):
			s, e = get_overlap((read.start1, read.stop1), (segment.reference_start, segment.reference_end));
			overlap1 = max(e-s, 0) 
		else:
			overlap1 = 0;
			
		if(read.chrom2 == segchrom and strand_conv[read.strand2] == segment.is_reverse):
			s, e = get_overlap((read.start2, read.stop2), (segment.reference_start, segment.reference_end));
			overlap2 = max(e-s, 0) 
		else:
			overlap2 = 0;
		
		return ComparisonChimera(rnum, read.stop1-read.start1, read.stop2-read.start2, read.typename, mapped_type, overlap1, overlap2);
			
	else:
		return ComparisonControl(rnum, read.length, read.typename, mapped_type, segment.get_tag('AS'))



samfile = pysam.Samfile(args.single_unique);
for segment in samfile.fetch(until_eof=True):
	rnum = int(segment.query_name.split("|")[0])
	segchrom = samfile.getrname(segment.tid)
	comparison_dict[rnum].append(single_unique_comparison(rnum, reads[rnum], segment, segchrom))
samfile.close();





#_____________________________________________________________________________________________________________________________
#control single reads
#_____________________________________________________________________________________________________________________________
#here overlap attribute of ComparisonChimera and ComparisonControl is an Alignment score of mapping
def single_control_comparison(rnum, read, segment):
	mapped_type = 'control'
	alignment_score = segment.get_tag('AS')
	if(read.typename == 'single'):
		return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type, alignment_score);
			
	elif(read.typename == 'chimera'):
		return ComparisonChimera(rnum, read.stop1-read.start1, read.stop2-read.start2, read.typename, mapped_type, alignment_score, alignment_score);
			
	else:
		return ComparisonControl(rnum, read.length, read.typename, mapped_type, alignment_score)
	
	

samfile = pysam.Samfile(args.single_control);
for segment in samfile.fetch(until_eof=True):
	rnum = int(segment.query_name.split("|")[0])
	comparison_dict[rnum].append(single_control_comparison(rnum, reads[rnum], segment))
samfile.close();






#_____________________________________________________________________________________________________________________________
#chimera unique reads
#_____________________________________________________________________________________________________________________________
def chimera_unique_comparison(rnum, read, segments, segchroms):
	mapped_type = 'chimera'
	if(read.typename == 'single'):
		for segment, segchrom in zip(segments, segchroms):
			if(read.chrom == segchrom and strand_conv[read.strand] == segment.is_reverse):
				s, e = get_overlap((read.start, read.stop), (segment.reference_start, segment.reference_end));
				overlap = max(e-s, 0)
				if(overlap>0):
					return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type, overlap);
		else:
			return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type, 0);
			
	elif(read.typename == 'chimera'):
		if(read.chrom1 == segchroms[0] and strand_conv[read.strand1] == segments[0].is_reverse):
			s, e = get_overlap((read.start1, read.stop1), (segments[0].reference_start, segments[0].reference_end));
			overlap1 = max(e-s, 0) 
		else:
			overlap1 = 0;
			
		if(read.chrom2 == segchroms[1] and strand_conv[read.strand2] == segments[1].is_reverse):
			s, e = get_overlap((read.start2, read.stop2), (segments[1].reference_start, segments[1].reference_end));
			overlap2 = max(e-s, 0) 
		else:
			overlap2 = 0;
		
		return ComparisonChimera(rnum, read.stop1-read.start1, read.stop2-read.start2, read.typename, mapped_type, overlap1, overlap2);
			
	else:
		return ComparisonControl(rnum, read.length, read.typename, mapped_type, max([x.get_tag('AS') for x in segments]))



samfile = pysam.Samfile(args.chimera_unique);
for segments in generator_doublesam(samfile):
	sort_segments_qstart(segments)
	rnum = int(segments[0].query_name.split("|")[0])
	segchroms = [samfile.getrname(segment.tid) for segment in segments]
	comparison_dict[rnum].append(chimera_unique_comparison(rnum, reads[rnum], segments, segchroms))
samfile.close();




#_____________________________________________________________________________________________________________________________
#control chimera reads
#_____________________________________________________________________________________________________________________________
#here overlap attribute of ComparisonChimera and ComparisonControl is a maximum(of two) Alignment score of mapping
def chimera_control_comparison(rnum, read, segments, segchroms):
	mapped_type_complete = 'control_chimera_complete'
	mapped_type_partial = 'control_chimera_partial'
	if(read.typename == 'single'):
		for segment, segchrom in zip(segments, segchroms):
			if(read.chrom == segchrom and strand_conv[read.strand] == segment.is_reverse):
				s, e = get_overlap((read.start, read.stop), (segment.reference_start, segment.reference_end));
				overlap = max(e-s, 0)
				if(overlap>0):
					return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type_partial, overlap);
		else:
			return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type_complete, 0);
			
	elif(read.typename == 'chimera'):
		if(read.chrom1 == segchroms[0] and strand_conv[read.strand1] == segments[0].is_reverse):
			s, e = get_overlap((read.start1, read.stop1), (segments[0].reference_start, segments[0].reference_end));
			overlap1 = max(e-s, 0) 
		else:
			overlap1 = 0;
			
		if(read.chrom2 == segchroms[1] and strand_conv[read.strand2] == segments[1].is_reverse):
			s, e = get_overlap((read.start2, read.stop2), (segments[1].reference_start, segments[1].reference_end));
			overlap2 = max(e-s, 0) 
		else:
			overlap2 = 0;
			
		if(overlap1>0 or overlap2>0):
			return ComparisonChimera(rnum, read.stop1-read.start1, read.stop2-read.start2, read.typename, mapped_type_partial, overlap1, overlap2);
		else:			
			return ComparisonChimera(rnum, read.stop1-read.start1, read.stop2-read.start2, read.typename, mapped_type_complete, overlap1, overlap2);			
			
	else:
		return ComparisonControl(rnum, read.length, read.typename, mapped_type_complete, max([x.get_tag('AS') for x in segments]))



samfile = pysam.Samfile(args.chimera_control);
for segments in generator_doublesam(samfile):
	sort_segments_qstart(segments)
	rnum = int(segments[0].query_name.split("|")[0])
	segchroms = [samfile.getrname(segment.tid) for segment in segments]
	comparison_dict[rnum].append(chimera_control_comparison(rnum, reads[rnum], segments, segchroms))
samfile.close();





#_____________________________________________________________________________________________________________________________
#nonunique single reads	
#_____________________________________________________________________________________________________________________________
def single_nonunique_comparison(rnum, read, segments, segchroms):
	mapped_type = 'single_nonunique'
	
	if(read.typename == 'single'):
		comparisons = [];
		for segment, segchrom in zip(segments, segchroms):
			comparisons.append(single_unique_comparison(rnum, read, segment, segchrom));
		best = max(comparisons, key=lambda x: x.overlap);
		return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type, best.overlap)

			
	elif(read.typename == 'chimera'):
		comparisons=[];
		for segment, segchrom in zip(segments, segchroms):
			comparisons.append(single_unique_comparison(rnum, read, segment, segchrom));
		best = max(comparisons, key=lambda x: x.overlap1+x.overlap2)	
		return ComparisonChimera(rnum, read.stop1-read.start1, read.stop2-read.start2, read.typename, mapped_type, best.overlap1, best.overlap2);

	else:
		return ComparisonControl(rnum, read.length, read.typename, mapped_type, max([x.get_tag('AS') for x in segments]))



samfile = pysam.Samfile(args.single_nonunique);
segments = [];
cnum = -1
for segment in samfile.fetch(until_eof=True):
	rnum = int(segment.query_name.split("|")[0])
	if(cnum==rnum):
		segments.append(segment);
	elif(segments):
		segchroms = [samfile.getrname(x.tid) for x in segments]
		comparison_dict[cnum].append(single_nonunique_comparison(cnum, reads[cnum], segments, segchroms))
		segments[:] = [segment]
		cnum = rnum
	else:
		cnum=rnum;
		segments.append(segment)
else:
	segchroms = [samfile.getrname(segment.tid) for segment in segments]
	comparison_dict[cnum].append(single_nonunique_comparison(cnum, reads[cnum], segments, segchroms))
samfile.close();





#_____________________________________________________________________________________________________________________________
#nonunique chimera reads	
#_____________________________________________________________________________________________________________________________
def chimera_nonunique_comparison(rnum, read, segment_blocks, segchrom_blocks):
	mapped_type = 'chimera_nonunique'
	
	if(read.typename == 'single'):
		comparisons=[];
		for segments, segchroms in zip(segment_blocks, segchrom_blocks):
			comparisons.append(chimera_unique_comparison(rnum, read, segments, segchroms));
		best = max(comparisons, key=lambda x: x.overlap)
		return ComparisonSingle(rnum, read.stop-read.start, read.typename, mapped_type, best.overlap)

	elif(read.typename == 'chimera'):
		comparisons=[];
		for segments, segchroms in zip(segment_blocks, segchrom_blocks):
			comparisons.append(chimera_unique_comparison(rnum, read, segments, segchroms));
		best = max(comparisons, key=lambda x: x.overlap1+x.overlap2)	
		return ComparisonChimera(rnum, read.stop1-read.start1, read.stop2-read.start2, read.typename, mapped_type, best.overlap1, best.overlap2);

	else:
		return ComparisonControl(rnum, read.length, read.typename, mapped_type, max([max(x[0].get_tag('AS'), x[1].get_tag('AS')) for x in segment_blocks]))


		
segment_blocks = [];
cnum = -1
samfile = pysam.Samfile(args.chimera_nonunique);
for segments in generator_doublesam(samfile):
	sort_segments_qstart(segments);
	rnum = int(segments[0].query_name.split("|")[0])
	if(cnum==rnum):
		segment_blocks.append(segments);
	elif(segment_blocks):
		segchrom_blocks = [(samfile.getrname(x[0].tid), samfile.getrname(x[1].tid)) for x in segment_blocks]
		comparison_dict[cnum].append(chimera_nonunique_comparison(cnum, reads[cnum], segment_blocks, segchrom_blocks))
		segment_blocks[:] = [segments]
		cnum = rnum
	else:
		cnum=rnum;
		segment_blocks.append(segments)
else:
	segchrom_blocks = [(samfile.getrname(x[0].tid), samfile.getrname(x[1].tid)) for x in segment_blocks]
	comparison_dict[cnum].append(chimera_nonunique_comparison(cnum, reads[cnum], segment_blocks, segchrom_blocks))	
samfile.close();		
			
			
			

#_____________________________________________________________________________________________________________________________			
#Get overall mapping stastics
def collapse_comparisons(comparisons):
	_order = {'single': 1, 'chimera': 2, 'control': 3, 'control_chimera_complete': 4, 'control_chimera_partial': 5, 'chimera_nonunique': 6, 'single_nonunique': 7}
	lc = len(comparisons);
	if(lc==1):
		return comparisons[0];
	elif(lc>1):	
		return min(comparisons, key=lambda x: _order[x.mapped_type])
		
			
def get_chimera_mapping(read, comparisons, fcutoff):
	if(comparisons):
		comp = collapse_comparisons(comparisons);
		is_correct = (float(comp.overlap1)/comp.length1 > fcutoff) and (float(comp.overlap2)/comp.length2 > fcutoff);
		return (comp.read_type, comp.mapped_type, is_correct, "|".join((read.ttype1, read.ttype2))), (comp.mapped_type, is_correct, comp.length1, comp.length2, read.conv1, read.conv2)
	else:
		return (read.typename, 'unmapped', False, "|".join((read.ttype1, read.ttype2))), ('unmapped', False, read.stop1-read.start1, read.stop2-read.start2, read.conv1, read.conv2)
	
	
def get_single_mapping(read, comparisons, fcutoff):
	if(comparisons):
		comp = collapse_comparisons(comparisons);
		is_correct = float(comp.overlap)/comp.length > fcutoff		
		return (comp.read_type, comp.mapped_type, is_correct, read.ttype), (comp.mapped_type, is_correct, comp.length, read.conv)
	else:
		return (read.typename, 'unmapped', False, read.ttype), ('unmapped', False, read.stop-read.start, read.conv)
		
		
def get_control_mapping(read, comparisons):
	if(comparisons):
		comp = collapse_comparisons(comparisons);
		return (comp.read_type, comp.mapped_type, False, read.ttype), (comp.mapped_type, False, comp.AS);
	else:
		return (read.typename, 'unmapped', False, read.ttype), ('unmapped', False, 0);
	
	

mapping_stat = defaultdict(int)
chimera_stat = defaultdict(int)
single_stat = defaultdict(int)
control_stat = defaultdict(int)
for num, read in enumerate(reads):
	if(read.typename=='chimera'):
		ms_key, ch_key = get_chimera_mapping(read, comparison_dict[num], args.coverage_cutoff)
		mapping_stat[ms_key] += 1;
		chimera_stat[ch_key] += 1
	elif(read.typename=='single'):
		ms_key, si_key = get_single_mapping(read, comparison_dict[num], args.coverage_cutoff)		
		mapping_stat[ms_key] += 1;
		single_stat[si_key] += 1;
	else:
		ms_key, co_key = get_control_mapping(read, comparison_dict[num])
		mapping_stat[ms_key] += 1;
		control_stat[co_key] += 1;
		
		
#for k, v in sorted(mapping_stat.items(), key=lambda x: x[1], reverse=True):
	#print "%s:\t%d" % (str(k), v)
	
	
if(not os.path.exists(args.outdir)):
    os.makedirs(args.outdir)	
	
with open(os.path.join(args.outdir, 'mapping_stat.yml'), 'w') as f:
	f.write(yaml.dump(dict(mapping_stat), default_flow_style=False))
	
with open(os.path.join(args.outdir, 'chimera_stat.yml'), 'w') as f:
	f.write(yaml.dump(dict(chimera_stat), default_flow_style=False))
	
with open(os.path.join(args.outdir, 'single_stat.yml'), 'w') as f:
	f.write(yaml.dump(dict(single_stat), default_flow_style=False))
	
with open(os.path.join(args.outdir, 'control_stat.yml'), 'w') as f:
	f.write(yaml.dump(dict(control_stat), default_flow_style=False))

	
#for num, comparisons in comparison_dict.items():
	#if(comparisons):
		#mapping_stat.append(collapse_comparisons(comparisons).mapped_type);
		##mapping_stat.extend([x.mapped_type for x in comparisons]);
		##print "\n".join([str(x) for x in comparisons])
		##print "_"*120
			
#print Counter(mapping_stat)
#print Counter([len(x) for x in comparison_dict.values()])
	
	
	

	
	
	