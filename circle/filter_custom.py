#! /usr/bin/python
'''Script tries to find reliable circles by applying following filters to the chimeric reads:
1) No more than 100 kilobase distance between chimeric hits
2) At least 2 independent reads reads supporing splice junction
3) GT/AG signal flanking the splice sites
4) Unambiguous chimeric breakpoint detection
''' 
import argparse
import sys;
from collections import defaultdict, Counter

from Bio import SeqIO
#from pybedtools import BedTool

from nrlbio.generators import generator_doublebed
from nrlbio.numerictools import distance as find_distance
from nrlbio.genome_system import seqrecord2seq
from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Script tries to find reliable circles by applying following filters to the chimeric reads:\n1) No more than 100 kilobase distance between chimeric hits\n2) At least 2 independent reads reads supporing splice junction\n3) GT/AG signal flanking the splice sites\n4) Unambiguous chimeric breakpoint detection');

parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the chimeras, double bed/gff file");
parser.add_argument('--distance', nargs = '?', default = 100000, type = int, help = "max allowed disctance between chimeric hits");
parser.add_argument('--reference', nargs = '?', required = True, type = str, help = "path to the reference(genome) to extract sequences from");
args = parser.parse_args();

#ASCONVERSION = {"+": "-", "-": "+"}

def check_distance(i1, i2, distance):
	return i1.chrom==i2.chrom and i1.strand==i2.strand and find_distance((i1.start, i1.end), (i2.start, i2.end))<distance;

def check_splice_site(i1, i2, seqrecord):
	gap = int(i1.attrs['gap'])
	if(i1.strand == '+'):
		flank5 = seqrecord2seq(seqrecord, i1.end+gap, i1.end+2, strand='+')
		flank3 = seqrecord2seq(seqrecord, i2.start-2, i2.start-gap, strand='+') 
		#if len(flank3)==2:
			#print flank3;
	elif(i1.strand == '-'):
		flank5 = seqrecord2seq(seqrecord, i1.start-2, i1.start-gap, strand='-') 
		flank3 = seqrecord2seq(seqrecord, i2.end+gap, i2.end+2, strand='-')
		#if len(flank3)==2:
			#print flank3;
	else:
		sys.stderr.write("Strand is not defined assumed to be plus\n");
		
	breakpoints = [];
	antisense = False;
	for b in range(abs(gap)+1):
		if(flank5[b:b+2] == "GT" and flank3[b:b+2] == "AG"):
			breakpoints.append(b);
		if(flank5[b:b+2] == "CT" and flank3[b:b+2] == "AC"):
			breakpoints.append(b);
			antisense = True;
	return breakpoints, antisense;
		
		
def clarify(i1, i2, breakpoint, antisense):
	gap = int(i1.attrs['gap']);
	#clarify breakpoint
	#clarify strandness
	if(antisense):
		if(i1.strand=='+'):
			cs1 = construct_gff_interval(chrom=i2.chrom, start=i2.start+breakpoint, stop=i2.stop, feature='ch', score=i2.score, strand='-', source='.', frame='.', attrs=i2.attrs.items())
			cs2 = construct_gff_interval(chrom=i1.chrom, start=i1.start, stop=i1.stop+gap+breakpoint, feature='ch', score=i1.score, strand='-', source='.', frame='.', attrs=i1.attrs.items())
		else:
			cs1 = construct_gff_interval(chrom=i2.chrom, start=i2.start, stop=i2.stop-breakpoint, feature='ch', score=i2.score, strand='+', source='.', frame='.', attrs=i2.attrs.items())
			cs2 = construct_gff_interval(chrom=i1.chrom, start=i1.start-gap-breakpoint, stop=i1.stop, feature='ch', score=i1.score, strand='+', source='.', frame='.', attrs=i1.attrs.items())
	else:
		if(i1.strand=='+'):
			cs1 = construct_gff_interval(chrom=i1.chrom, start=i1.start, stop=i1.stop+gap+breakpoint, feature='ch', score=i1.score, strand=i1.strand, source='.', frame='.', attrs=i1.attrs.items())
			cs2 = construct_gff_interval(chrom=i2.chrom, start=i2.start+breakpoint, stop=i2.stop, feature='ch', score=i2.score, strand=i2.strand, source='.', frame='.', attrs=i2.attrs.items())
		else:	
			cs1 = construct_gff_interval(chrom=i1.chrom, start=i1.start-gap-breakpoint, stop=i1.stop, feature='ch', score=i1.score, strand=i1.strand, source='.', frame='.', attrs=i1.attrs.items())
			cs2 = construct_gff_interval(chrom=i2.chrom, start=i2.start, stop=i2.stop-breakpoint, feature='ch', score=i2.score, strand=i2.strand, source='.', frame='.', attrs=i2.attrs.items())
		
	return cs1, cs2;


def assign_splice_type(i1, i2):
	if(i1.strand=='+'):
		if(i1.start<i2.start):
			return 'lsj';
		else:
			return 'csj';
	else:
		if(i1.start>i2.start):
			return 'lsj';
		else:
			return 'csj';

		
	
	
	
reference = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

total = 0;
passed = 0;
uniqbreakpoints = [];
stypes = defaultdict(int);
for i1, i2 in generator_doublebed(args.path):
	total += 1;
	if(check_distance(i1, i2, args.distance)):
		seqrecord = reference[i1.chrom];
		breakpoints, antisense = check_splice_site(i1, i2, seqrecord);
		uniqbreakpoints.append(len(breakpoints))
		if(len(breakpoints)==1):
			cs1, cs2 = clarify(i1, i2, breakpoints[0], antisense)
			splice_type = assign_splice_type(cs1, cs2)
			passed += 1;
		else:
			splice_type = 'intra'
			cs1, cs2 = i1, i2
	else:
		splice_type = 'inter'
		cs1, cs2 = i1, i2
		
	cs1.attrs['interaction'] = splice_type
	cs2.attrs['interaction'] = splice_type
	sys.stdout.write("%s%s" % (cs1, cs2));	
	stypes[splice_type]+=1;
		
		

num_splice_sites = Counter(uniqbreakpoints);
sys.stderr.write("total chimeras: %d\npassed chimeras: %d\nfraction passed %1.5f\n\n" % (total, passed, float(passed)/total));
sys.stderr.write("Ambiguity	in breackpoint detection\nnum canonical splice sites\tnum chimeras\n");
for k in sorted(num_splice_sites.keys()):
	sys.stderr.write("%d\t%d\n" % (k, num_splice_sites[k]));
sys.stderr.write("\nSplice site type\tnumber of reads\n");
for k, v in stypes.items():
	sys.stderr.write("%s\t%d\n" % (k, v));






