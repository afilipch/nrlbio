#! /usr/lib/python
'''Generates "real" and "fake" chimeras for the sake of testing''' 
import argparse
import sys;
import random;
from collections import defaultdict, namedtuple

from pybedtools import Interval
from Bio import SeqIO

from nrlbio.genome_system import generate_genes_from_refseq;
from nrlbio.sequencetools import random_string, shuffle_string
from nrlbio.random_extension import weighted2interval, weighted_choice_fast


parser = argparse.ArgumentParser(description='Generates "real" and "fake" chimeras for the sake of testing');
#parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file to be annotated");
parser.add_argument('-f', '--fasta', nargs = '?', required = True, type = str, help = "path to a fasta file to generate \"real\" chimeras from");
parser.add_argument('-g', '--genes', nargs = '?', required = True, type = str, help = "path to a genome system file(refseq genes)");
parser.add_argument('-l', '--length', nargs = '?', default = 100, type = int, help = "length of generated chimeras");
parser.add_argument('-mfl', '--min_fragment_length', nargs = '?', default = 1, type = int, help = "min length of chimera parts");
parser.add_argument('-ncr', '--num_chimeric_reads', nargs = 6, default = (10000, 5000, 5000, 10000, 5000, 10000), type = int, help = "number of exonic-exonic, exonic-intronic, exonic-intergenic, intronic-intronic, intronic-intergenic, intergenic-intergenic chimeras");
parser.add_argument('-nsr', '--num_single_reads', nargs = 3, default = (100000, 50000, 100000), type = int, help = "number of exonic, intronic, intergenic single reads");
parser.add_argument('-nrr', '--num_random_reads', nargs = 2, default = (20000, 60000), type = int, help = "number of shuffled, random reads");
parser.add_argument('-lr', '--left_random', nargs = 2, default = (1, 0), type = float, help = "options of 5' random sequence length (Probability of 0, maximum length of the random head)");
parser.add_argument('-rr', '--right_random', nargs = 2, default = (1, 0), type = float, help = "options of 3' random sequence length (Probability of 0, maximum length of the random tail)");
args = parser.parse_args();

# set left and right random sequence length probability
def _get_intervals_items(zero_prob, max_length):
	max_length = int(max_length)
	l = [(0, zero_prob)];
	if(max_length):
		uniform_prob = (1-zero_prob)/max_length
		for n in range(1, max_length+1):
			l.append((n, uniform_prob));
	return	weighted2interval(l)

litems, lintervals = _get_intervals_items(*args.left_random)
ritems, rintervals = _get_intervals_items(*args.right_random)
#_______________________________________________________________________________________________________________



def select_exon(genes, length):
	while(True):
		gene = random.choice(genes)
		exon = random.choice(gene.exons)
		if(len(exon)>=length):
			return exon;

def select_intron(genes, length):
	while(True):
		gene = random.choice(genes)
		if(gene.introns):
			intron = random.choice(gene.introns)
			if(len(intron)>=length):
				return intron;

def select_chromosome(chromosomes, length):
	while(True):
		chrom = random.choice(chromosomes)
		if(len(chrom)>=length):
			return chrom
#_______________________________________________________________________________________________________________



def get_random_interval(interval, length):
	start = random.randint(0, len(interval) - length) + interval.start
	end = start + length;
	return Interval(interval.chrom, start, end, name=".", score=".", strand=interval.strand, otherfields=None)


def get_random_intervals(i1, i2, length, joint = 0, mfl = 1):
	if(not joint):
		joint = random.randint(mfl, length-(mfl+1));
	return get_random_interval(i1, joint), get_random_interval(i2, (length - joint))	


def get_seq(interval, fasta_dict, left, right):
	if(interval.strand == "+"):
		seq = str(fasta_dict[interval.chrom].seq[interval.start:interval.end].upper());
	else:
		seq = str(fasta_dict[interval.chrom].seq[interval.start:interval.end].reverse_complement().upper());	
		
	if("N" in seq):
		return None;
		
	lseq = random_string(left);
	rseq = random_string(right);
	
	return "".join((lseq, seq, rseq));
#_______________________________________________________________________________________________________________

	
	
def interval2read(interval, fasta_dict, type_, left=0, right=0):
	header = "|".join(["single", '%s:%s:%d:%d' % (interval.chrom, interval.strand, interval.start, interval.end), type_, '%d:%d' % (left, right)])
	seq = get_seq(interval, fasta_dict, left, right)
	
	if(seq):
		return ">%s\n%s" % (header, seq);
	else:
		return None
	
	
def interval2shuffled(interval, fasta_dict, type_):
	header = "|".join(["shuffled", '%s:%s:%d:%d' % (interval.chrom, interval.strand, interval.start, interval.end), type_, str(len(interval))])
	seq = get_seq(interval, fasta_dict, 0, 0)
	
	if(seq):
		return ">%s\n%s" % (header, shuffle_string(seq));
	else:
		return None
	
	
def get_random_entry(length):
	header = "|".join(["random", '0', '0', str(length)])
	seq = random_string(length);	
	return ">%s\n%s" % (header, seq);

	
def intervals2chimera(intervals, fasta_dict, type_, gap = 0, left=0, right=0):
	header = "|".join(['chimera', ";".join(['%s:%s:%d:%d' % (x.chrom, x.strand, x.start, x.end) for x in intervals]), type_, '%d:%d:%d' % (left, right, gap)])
	seqs = get_seq(intervals[0], fasta_dict, left, gap), get_seq(intervals[1], fasta_dict, 0, right)
	
	if(all(seqs)):
		return ">%s\n%s" % (header, "".join(seqs))
	else:
		return None


	

	

# convert input instances into python objects
fasta = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))

chromosomes = [];
for k, v in fasta.iteritems():
	for strand in "+-":
		chromosomes.append(  Interval(k, 0, len(v), name=".", score=".", strand=strand, otherfields=None) )

genes = list(generate_genes_from_refseq(args.genes));
#for c, gene in enumerate(generate_genes_from_refseq(args.genes)):
	#genes[c, len(gene.transcript)] = gene;
#_______________________________________________________________________________________________________________


#exonic-exonic
for _ in range(args.num_chimeric_reads[0]):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		exon1 = select_exon(genes, joint)
		exon2 = select_exon(genes, args.length-joint)
		entry = intervals2chimera(get_random_intervals(exon1, exon2, args.length, joint = joint), fasta, "exon:exon")
		if(entry):
			print entry;
			break;
#_______________________________________________________________________________________________________________


#exonic-intronic
for _ in range(args.num_chimeric_reads[1]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		exon = select_exon(genes, joint)
		intron = select_intron(genes, args.length-joint)
		entry = intervals2chimera(get_random_intervals(exon, intron, args.length, joint = joint), fasta, "exon:intron")
		if(entry):
			print entry;
			break;

#_______________________________________________________________________________________________________________


#intronic-exonic
for _ in range(args.num_chimeric_reads[1]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		intron = select_intron(genes, joint)
		exon = select_exon(genes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(intron, exon, args.length, joint = joint), fasta, "intron:exon")
		if(entry):
			print entry;
			break;
		
#_______________________________________________________________________________________________________________


#exonic-intergenic
for _ in range(args.num_chimeric_reads[2]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		exon = select_exon(genes, joint)
		chrom = select_chromosome(chromosomes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(exon, chrom, args.length, joint = joint), fasta, "exon:intergenic")
		if(entry):
			print entry;
			break;
		
#_______________________________________________________________________________________________________________


#intergenic-exonic
for _ in range(args.num_chimeric_reads[2]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		chrom = select_chromosome(chromosomes, joint)
		exon = select_exon(genes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(chrom, exon, args.length, joint = joint), fasta, "intergenic:exon")
		if(entry):
			print entry;
			break;

#_______________________________________________________________________________________________________________


#intronic-intronic
for _ in range(args.num_chimeric_reads[3]):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		intron1 = select_intron(genes, joint)
		intron2 = select_intron(genes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(intron1, intron2, args.length, joint = joint), fasta, "intron:intron")
		if(entry):
			print entry;
			break;

#_______________________________________________________________________________________________________________


#intronic-intergenic
for _ in range(args.num_chimeric_reads[4]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		intron = select_intron(genes, joint)
		chrom = select_chromosome(chromosomes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(intron, chrom, args.length, joint = joint), fasta, "intron:intergenic")
		if(entry):
			print entry;
			break;
		

#_______________________________________________________________________________________________________________


#intergenic-intronic
for _ in range(args.num_chimeric_reads[4]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		chrom = select_chromosome(chromosomes, joint)
		intron = select_intron(genes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(chrom, intron, args.length, joint = joint), fasta, "intergenic:intron")
		if(entry):
			print entry;
			break;

#_______________________________________________________________________________________________________________


#intergenic-intergenic
for _ in range(args.num_chimeric_reads[5]):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
		
	while(True):
		chrom1 = select_chromosome(chromosomes, joint)
		chrom2 = select_chromosome(chromosomes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(chrom1, chrom2, args.length, joint = joint), fasta, "intergenic:intergenic")
		if(entry):
			print entry;
			break;
#_______________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________




#exonic single reads
for _ in range(args.num_single_reads[0]):
	left = weighted_choice_fast(litems, lintervals)
	right = weighted_choice_fast(ritems, rintervals)
	length = args.length - left - right
	while(True):
		exon = select_exon(genes, length)
		entry = interval2read(get_random_interval(exon, length), fasta, "exon", left=left, right=right)
		if(entry):
			print entry;
			break;
#_______________________________________________________________________________________________________________


#intronic single reads
for _ in range(args.num_single_reads[1]):
	left = weighted_choice_fast(litems, lintervals)
	right = weighted_choice_fast(ritems, rintervals)
	length = args.length - left - right
	while(True):
		exon = select_intron(genes, length)
		entry = interval2read(get_random_interval(intron, length), fasta, "intron", left=left, right=right)
		if(entry):
			print entry;
			break;
#_______________________________________________________________________________________________________________


#intergenic single reads
for _ in range(args.num_single_reads[2]):
	left = weighted_choice_fast(litems, lintervals)
	right = weighted_choice_fast(ritems, rintervals)
	length = args.length - left - right
	while(True):
		chrom = select_chromosome(chromosomes, length)
		entry = interval2read(get_random_interval(chrom, length), fasta, "intergenic", left=left, right=right)
		if(entry):
			print entry;
			break;
#_______________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________



#shuffled exonic reads
for _ in range(args.num_random_reads[0]/3):
	while(True):
		exon = select_exon(genes, args.length)		
		entry = interval2shuffled(get_random_interval(exon, args.length), fasta, "exon")
		if(entry):
			print entry;
			break;
#_______________________________________________________________________________________________________________


#shuffled intronic reads
for _ in range(args.num_random_reads[0]/3):
	while(True):
		exon = select_intron(genes, args.length)	
		entry = interval2shuffled(get_random_interval(intron, args.length), fasta, "intron")
		if(entry):
			print entry;
			break;
#_______________________________________________________________________________________________________________


#shuffled intergenic reads:
for _ in range(args.num_random_reads[0]/3):
	while(True):
		chrom = select_chromosome(chromosomes, args.length)
		entry = interval2shuffled(get_random_interval(chrom, args.length), fasta, "intergenic")
		if(entry):
			print entry;
			break;
#_______________________________________________________________________________________________________________

#random reads:
for _ in range(args.num_random_reads[1]):
	print get_random_entry(args.length);




























	
	
	