#! /usr/lib/python
'''Generates "real" and "fake" chimeras for the sake of testing''' 
import argparse
import sys;
import random;
from collections import defaultdict, namedtuple

from pybedtools import Interval
from Bio import SeqIO

from nrlbio.genome_system import generate_transcripts_from_refseq;
from nrlbio.sequencetools import random_string, shuffle_string, introduce_conversions
from nrlbio.random_extension import weighted2interval, weighted_choice_fast


parser = argparse.ArgumentParser(description='Generates "real" and "fake" chimeras for the sake of testing');
#parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file to be annotated");
parser.add_argument('-f', '--fasta', nargs = '?', required = True, type = str, help = "path to a fasta file to generate \"real\" chimeras from");
parser.add_argument('-g', '--genes', nargs = '?', required = True, type = str, help = "path to a genome system file(refseq genes)");

parser.add_argument('-l', '--length', nargs = '?', default = 100, type = int, help = "length of generated chimeras");
parser.add_argument('-mfl', '--min_fragment_length', nargs = '?', default = 1, type = int, help = "min length of chimera parts");

parser.add_argument('-ncr', '--num_chimeric_reads', nargs = 6, default = (5000, 5000, 5000, 5000, 5000, 5000), type = int, help = "number of exonic-exonic, exonic-intronic, exonic-intergenic, intronic-intronic, intronic-intergenic, intergenic-intergenic chimeras");
parser.add_argument('-nsr', '--num_single_reads', nargs = 3, default = (100000, 100000, 100000), type = int, help = "number of exonic, intronic, intergenic single reads");
parser.add_argument('-nrr', '--num_random_reads', nargs = 2, default = (150000, 150000), type = int, help = "number of shuffled, random reads");

parser.add_argument('-lr', '--left_random', nargs = 2, default = (1, 0), type = float, help = "options of 5' random sequence length (Probability of 0, maximum length of the random head)");
parser.add_argument('-rr', '--right_random', nargs = 2, default = (1, 0), type = float, help = "options of 3' random sequence length (Probability of 0, maximum length of the random tail)");

parser.add_argument('--conv_prob', nargs = '?', default = 0, type = float, help = "probability of conversion");

parser.add_argument('--mir', nargs = '?', type = str, help = "path to a fasta file with miRNAs. If not set miRNA chimeras will not be generated");
parser.add_argument('-nmr', '--num_mir_reads', nargs = 6, default = (50000, 5000, 5000, 5000, 5000, 5000), type = int, help = "number of miRNAs only, miRNA-exonic chimeras, miRNA-intronic chimeras, miRNA-intergenic chimeras, miRNA-shuffled chimeras, miRNA-random chimeras to be generated");
parser.add_argument('-mlr', '--mir_left_random', nargs = 2, default = (0.9, 4), type = float, help = "options of 5' random sequence length for miRNAs and miRNA chimeras (Probability of 0, maximum length of the random head)");
parser.add_argument('-mrr', '--mir_right_random', nargs = 2, default = (0.5, 10), type = float, help = "options of 3' random sequence length for miRNAs and miRNA chimeras (Probability of 0, maximum length of the random tail)");
parser.add_argument('-mlc', '--mir_left_cut', nargs = 2, default = (0.9, 4), type = float, help = "options of 5' cut length for miRNAs and miRNA chimeras (Probability of 0, maximum length of the random cut)");
parser.add_argument('-mrc', '--mir_right_cut', nargs = 2, default = (0.1, 8), type = float, help = "options of 3' cut length for miRNAs and miRNA chimeras (Probability of 0, maximum length of the random cut)");
parser.add_argument('--target_length', nargs = '?', default = 50, type = int, help = "length of mirna target sequences");

args = parser.parse_args();

TOTAL_COUNT=0#count of reads generated
ProxyInterval = namedtuple('Proxy', ['chrom', 'strand', 'start', 'end'])
RANDOM_INTERVAL = ProxyInterval('chr', '+', 10, 20)


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


def get_seq(interval, fasta_dict, left, right, conv_prob):
	if(interval.strand == "+"):
		seq = str(fasta_dict[interval.chrom].seq[interval.start:interval.end].upper());
	else:
		seq = str(fasta_dict[interval.chrom].seq[interval.start:interval.end].reverse_complement().upper());	
	
	#sequences with N(undefined nucleotis) are skipped for following reasons: 1) if N is inside the seq it would harm mapping performance. 2) Reads with N are anyway removed during real read preproccessing step
	if("N" in seq):
		return None, 0
	
	#add random flnaking sequences (unremoved adapters or sequencing machine additions)
	lseq = random_string(left);
	rseq = random_string(right);
	
	#introduce conversions
	if(conv_prob):
		seq, conv_num = introduce_conversions(seq, conv_prob)
	else:
		conv_num = 0;
		
	return "".join((lseq, seq, rseq)), conv_num;
#_______________________________________________________________________________________________________________

	
	
def interval2read(interval, fasta_dict, type_, left=0, right=0, conv_prob=0, type_name='single'):
	seq, conv_num = get_seq(interval, fasta_dict, left, right, conv_prob)
	header = "|".join([str(TOTAL_COUNT), type_name, '%s:%s:%d:%d' % (interval.chrom, interval.strand, interval.start, interval.end), type_, '%d:%d:%d' % (left, right, 0), str(conv_num)])
	
	if(seq):
		return ">%s\n%s" % (header, seq);
	else:
		return None
	
	
def interval2shuffled(interval, fasta_dict, type_):
	header = "|".join([str(TOTAL_COUNT), "shuffled", '%s:%s:%d:%d' % (interval.chrom, interval.strand, interval.start, interval.end), type_, '%d:%d:%d' % (len(interval), 0, 0), '0'])
	seq, stub = get_seq(interval, fasta_dict, 0, 0, 0)
	
	if(seq):
		return ">%s\n%s" % (header, shuffle_string(seq));
	else:
		return None
	
	
def get_random_entry(length):
	header = "|".join([str(TOTAL_COUNT), "random", '%s:%s:%d:%d' % (RANDOM_INTERVAL.chrom, RANDOM_INTERVAL.strand, RANDOM_INTERVAL.start, RANDOM_INTERVAL.end), 'random', '%d:%d:%d' % (length, 0, 0), '0'])
	seq = random_string(length);
	return ">%s\n%s" % (header, seq);

	
def intervals2chimera(intervals, fasta_dict, type_, gap = 0, left=0, right=0, conv_prob=0):
	seq1, conv_num1 = get_seq(intervals[0], fasta_dict, left, gap, conv_prob)
	seq2, conv_num2 = get_seq(intervals[1], fasta_dict, 0, right, conv_prob)
	header = "|".join([str(TOTAL_COUNT), 'chimera', "&".join(['%s:%s:%d:%d' % (x.chrom, x.strand, x.start, x.end) for x in intervals]), type_, '%d:%d:%d' % (left, right, gap), "%d:%d" % (conv_num1, conv_num2)])
	
	if(seq1 and seq2):
		return ">%s\n%s" % (header, "".join((seq1, seq2)))
	else:
		return None


	

	

# convert input instances into python objects
fasta = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))

chromosomes = [];
for k, v in fasta.iteritems():
	for strand in "+-":
		chromosomes.append(  Interval(k, 0, len(v), name=".", score=".", strand=strand, otherfields=None) )

genes = list(generate_transcripts_from_refseq(args.genes));
#for c, gene in enumerate(generate_transcripts_from_refseq(args.genes)):
	#genes[c, len(gene.transcript)] = gene;
	
	

#_______________________________________________________________________________________________________________


#exonic-exonic
for _ in range(args.num_chimeric_reads[0]):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		exon1 = select_exon(genes, joint)
		exon2 = select_exon(genes, args.length-joint)
		entry = intervals2chimera(get_random_intervals(exon1, exon2, args.length, joint = joint), fasta, "exon:exon", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;
#_______________________________________________________________________________________________________________


#exonic-intronic
for _ in range(args.num_chimeric_reads[1]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		exon = select_exon(genes, joint)
		intron = select_intron(genes, args.length-joint)
		entry = intervals2chimera(get_random_intervals(exon, intron, args.length, joint = joint), fasta, "exon:intron", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;

#_______________________________________________________________________________________________________________


#intronic-exonic
for _ in range(args.num_chimeric_reads[1]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		intron = select_intron(genes, joint)
		exon = select_exon(genes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(intron, exon, args.length, joint = joint), fasta, "intron:exon", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;
		
#_______________________________________________________________________________________________________________


#exonic-intergenic
for _ in range(args.num_chimeric_reads[2]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		exon = select_exon(genes, joint)
		chrom = select_chromosome(chromosomes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(exon, chrom, args.length, joint = joint), fasta, "exon:intergenic", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;
		
#_______________________________________________________________________________________________________________


#intergenic-exonic
for _ in range(args.num_chimeric_reads[2]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		chrom = select_chromosome(chromosomes, joint)
		exon = select_exon(genes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(chrom, exon, args.length, joint = joint), fasta, "intergenic:exon", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;

#_______________________________________________________________________________________________________________


#intronic-intronic
for _ in range(args.num_chimeric_reads[3]):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		intron1 = select_intron(genes, joint)
		intron2 = select_intron(genes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(intron1, intron2, args.length, joint = joint), fasta, "intron:intron", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;

#_______________________________________________________________________________________________________________


#intronic-intergenic
for _ in range(args.num_chimeric_reads[4]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		intron = select_intron(genes, joint)
		chrom = select_chromosome(chromosomes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(intron, chrom, args.length, joint = joint), fasta, "intron:intergenic", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;
		

#_______________________________________________________________________________________________________________


#intergenic-intronic
for _ in range(args.num_chimeric_reads[4]/2):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
	
	while(True):
		chrom = select_chromosome(chromosomes, joint)
		intron = select_intron(genes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(chrom, intron, args.length, joint = joint), fasta, "intergenic:intron", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;

#_______________________________________________________________________________________________________________


#intergenic-intergenic
for _ in range(args.num_chimeric_reads[5]):
	joint = random.randint(args.min_fragment_length, args.length-(args.min_fragment_length +1));
		
	while(True):
		chrom1 = select_chromosome(chromosomes, joint)
		chrom2 = select_chromosome(chromosomes, args.length - joint)
		entry = intervals2chimera(get_random_intervals(chrom1, chrom2, args.length, joint = joint), fasta, "intergenic:intergenic", conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
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
		entry = interval2read(get_random_interval(exon, length), fasta, "exon", left=left, right=right, conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;
#_______________________________________________________________________________________________________________


#intronic single reads
for _ in range(args.num_single_reads[1]):
	left = weighted_choice_fast(litems, lintervals)
	right = weighted_choice_fast(ritems, rintervals)
	length = args.length - left - right
	while(True):
		intron = select_intron(genes, length)
		entry = interval2read(get_random_interval(intron, length), fasta, "intron", left=left, right=right, conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
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
		entry = interval2read(get_random_interval(chrom, length), fasta, "intergenic", left=left, right=right, conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
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
			TOTAL_COUNT+=1;
			print entry;
			break;
#_______________________________________________________________________________________________________________


#shuffled intronic reads
for _ in range(args.num_random_reads[0]/3):
	while(True):
		exon = select_intron(genes, args.length)
		entry = interval2shuffled(get_random_interval(intron, args.length), fasta, "intron")
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;
#_______________________________________________________________________________________________________________


#shuffled intergenic reads:
for _ in range(args.num_random_reads[0]/3):
	while(True):
		chrom = select_chromosome(chromosomes, args.length)
		entry = interval2shuffled(get_random_interval(chrom, args.length), fasta, "intergenic")
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;
#_______________________________________________________________________________________________________________

#random reads:
for _ in range(args.num_random_reads[1]):
	print get_random_entry(args.length);
	TOTAL_COUNT+=1;







#_______________________________________________________________________________________________________________
#miRNA chimeras part
#_______________________________________________________________________________________________________________
#get mirna dictionary from fasta file
if(args.mir):
	mirdict = SeqIO.to_dict(SeqIO.parse(args.mir, 'fasta'));
else:
	sys.exit();

# set left and right random sequence length probability
litems_cut, lintervals_cut = _get_intervals_items(*args.mir_left_cut)
ritems_cut, rintervals_cut = _get_intervals_items(*args.mir_right_cut)
litems_mir, lintervals_mir = _get_intervals_items(*args.mir_left_random)
ritems_mir, rintervals_mir = _get_intervals_items(*args.mir_right_random)

def _get_mir_edges():
	return weighted_choice_fast(litems_mir, lintervals_mir), weighted_choice_fast(ritems_mir, rintervals_mir), weighted_choice_fast(litems_cut, lintervals_cut), weighted_choice_fast(ritems_cut, rintervals_cut)

def _get_mir_interval(mirname, mirseq, left_cut, right_cut):
	return Interval(mirname, left_cut, len(mirseq)-right_cut, name=".", score=".", strand='+', otherfields=None)


def get_mir_chimera(m_int, t_int, mirdict, target_dict, type_, gap = 0, left=0, right=0, conv_prob=0):
	
	seq1, conv_num1 = get_seq(m_int, mirdict, left, gap, conv_prob)
	seq2, conv_num2 = get_seq(t_int, target_dict, 0, right, conv_prob)
	header = "|".join([str(TOTAL_COUNT), 'mirna_chimera', "&".join(['%s:%s:%d:%d' % (x.chrom, x.strand, x.start, x.end) for x in (m_int, t_int)]), type_, '%d:%d:%d' % (left, right, gap), "%d:%d" % (conv_num1, conv_num2)])
	
	if(seq1 and seq2):
		return ">%s\n%s" % (header, "".join((seq1, seq2)))
	else:
		return None
	
	
def get_mir_shuffled(m_int, t_int, mirdict, target_dict, type_, left=0, conv_prob=0):	
	
	seq1, conv_num1 = get_seq(m_int, mirdict, left, 0, conv_prob)
	seq2 = shuffle_string(get_seq(t_int, target_dict, 0, 0, 0)[0])
	header = "|".join([str(TOTAL_COUNT), 'mirna_shuffled', "&".join(['%s:%s:%d:%d' % (x.chrom, x.strand, x.start, x.end) for x in (m_int, t_int)]), type_, '%d:%d:%d' % (left, 0, 0), "%d:%d" % (conv_num1, 0)])
	
	if(seq1 and seq2):
		return ">%s\n%s" % (header, "".join((seq1, seq2)))
	else:
		return None
	
	
def get_mir_random(m_int, mirdict, length, type_, left=0, conv_prob=0):	
	
	seq1, conv_num1 = get_seq(m_int, mirdict, left, 0, conv_prob)
	seq2 = random_string(length);
	header = "|".join([str(TOTAL_COUNT), 'mirna_random', "&".join(['%s:%s:%d:%d' % (x.chrom, x.strand, x.start, x.end) for x in (m_int, RANDOM_INTERVAL)]), type_, '%d:%d:%d' % (left, length,0), "%d:%d" % (conv_num1, 0)])
	
	if(seq1 and seq2):
		return ">%s\n%s" % (header, "".join((seq1, seq2)))
	else:
		return None
	
	
	
	
	

#miRNA single reads
for _ in range(args.num_mir_reads[0]):
	left, right, left_cut, right_cut = _get_mir_edges();
	
	mirname = random.choice(mirdict.keys())
	mirseq = mirdict[mirname]	
	mir_interval = _get_mir_interval(mirname, mirseq, left_cut, right_cut)

	print interval2read(mir_interval, mirdict, "mirna", left=left, right=right, conv_prob=args.conv_prob, type_name='mirna_single')
	TOTAL_COUNT+=1;




#miRNA chimera exon
for _ in range(args.num_mir_reads[1]):
	left, right, left_cut, right_cut = _get_mir_edges();
	right_target = weighted_choice_fast(ritems, rintervals)
	
	mirname = random.choice(mirdict.keys())
	mirseq = mirdict[mirname]	
	mir_interval = _get_mir_interval(mirname, mirseq, left_cut, right_cut)
	
	while(True):
		exon = select_exon(genes, args.target_length)
		entry = get_mir_chimera(mir_interval, get_random_interval(exon, args.target_length), mirdict, fasta, "exon", gap = right, left=left, right=right_target, conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;



#miRNA chimera intron
for _ in range(args.num_mir_reads[2]):
	left, right, left_cut, right_cut = _get_mir_edges();
	right_target = weighted_choice_fast(ritems, rintervals)
	
	mirname = random.choice(mirdict.keys())
	mirseq = mirdict[mirname]	
	mir_interval = _get_mir_interval(mirname, mirseq, left_cut, right_cut)
	
	while(True):
		intron = select_intron(genes, args.target_length)
		entry = get_mir_chimera(mir_interval, get_random_interval(intron, args.target_length), mirdict, fasta, "intron", gap = right, left=left, right=right_target, conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;



#miRNA chimera intergenic
for _ in range(args.num_mir_reads[3]):
	left, right, left_cut, right_cut = _get_mir_edges();
	right_target = weighted_choice_fast(ritems, rintervals)
	
	mirname = random.choice(mirdict.keys())
	mirseq = mirdict[mirname]	
	mir_interval = _get_mir_interval(mirname, mirseq, left_cut, right_cut)
	
	while(True):
		chrom = select_chromosome(chromosomes, args.target_length)
		entry = get_mir_chimera(mir_interval, get_random_interval(chrom, args.target_length), mirdict, fasta, "intergenic", gap = right, left=left, right=right_target, conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;



#miRNA chimera exon shuffled
for _ in range(args.num_mir_reads[4]):
	left, right, left_cut, right_cut = _get_mir_edges();
	
	mirname = random.choice(mirdict.keys())
	mirseq = mirdict[mirname]	
	mir_interval = _get_mir_interval(mirname, mirseq, left_cut, right_cut)
	
	while(True):
		exon = select_exon(genes, args.target_length)
		entry = get_mir_shuffled(mir_interval, get_random_interval(exon, args.target_length), mirdict, fasta, "exon", left=left, conv_prob=args.conv_prob)
		if(entry):
			TOTAL_COUNT+=1;
			print entry;
			break;


#miRNA chimera exon random
for _ in range(args.num_mir_reads[5]):
	left, right, left_cut, right_cut = _get_mir_edges();

	mirname = random.choice(mirdict.keys())
	mirseq = mirdict[mirname]	
	mir_interval = _get_mir_interval(mirname, mirseq, left_cut, right_cut)
	
	print get_mir_random(mir_interval, mirdict, args.target_length, "random", left=left, conv_prob=args.conv_prob)
	TOTAL_COUNT+=1;


	
	