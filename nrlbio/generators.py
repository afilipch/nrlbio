# /usr/bin/python
'''generators for different kinds of input files. They may be used to allow multiprocessing and save memory footprint.
usually you should provide path to input file and size of elements(not lines) in buffer'''
import sys;
import re;
import itertools;
import sequencetools
from collections import *


fastq = namedtuple('fastq', 'id, seq, sign, qual')

def generator_fastq(path, take = ['id', 'seq', 'sign', 'qual'], reverse = False, shuffle = False):
	'''yields fastq object(collections.namedtuple based) from fastq file provided
	
	path str: path to fastq file
	take list: list of fastq entry attributes to pass to yielded fastq object. Allows to save some memory. Entries have to be: "id", "seq", "sign", "qual"
	reverse bool: if True fastq sequence and quality string will be reversed
	shuffle bool: if True fastq sequence will be reversed
	
	Yields fastq: fastq object(collections.namedtuple based) 
	'''
	my_take = {"id": 0, "seq": 1, "sign": 2, "qual": 3}
	line_set = [];
	
	with open(path) as f:
		for line in f:  
			line_set.append(line.strip());
			if (len(line_set) == 4):
				appendix = [None]*4;
				for el in take:
					if(reverse and el in ["seq", "qual"]):
						appendix[my_take[el]] = line_set[my_take[el]][::-1];
					elif(shuffle and el == "seq"):
							appendix[my_take[el]] = sequencetools.shuffle_string(line_set[my_take[el]])
					else:	
						appendix[my_take[el]] = line_set[my_take[el]]
				yield fastq(*appendix)
				line_set = [];				

	
def generator_maf(path, aligned_species=None):
	'''reads maf stitched file, taking into account interactions id. Yields conservation.Maf object corresponding to each entry in maf file
	
	path string: path to maf stitched file
	aligned_species list: names of aligned species(genome suystems: ce6, mm9, etc.) to be considered in further analysis. If not given takes all species
	
	Yields conservation.Maf: object corresponding to each entry in maf file
	'''	
	
	from conservation import Maf;
	
	nonaligned = re.compile('-+$')
	alignment = OrderedDict()
	with open(path) as f:
		for l in f:
			l = l.strip();  
			if(l):
				if(l[0] == ">"):
					if(len(l.split(":")) == 3):
						header = l
						mflag = True;
					else:
						specie = l[1:];
						mflag = False;
				elif(aligned_species == None or specie in aligned_species):
					if(mflag):
						refseq = l.upper();
					elif(not nonaligned.match(l)):
						alignment[specie] = l.upper();
			else:
				yield Maf(header, alignment, refseq)
		try:
			yield Maf(header, alignment)
		except:
			pass
		
		
def generator_doublebed(path):
	'''Reads \'double\' bed file. Yields consecutive pairs of bed/gff intervals
	
	path string: path to maf stitched file
	
	Yields list: 2-element list of BedTool intervals, representing chimera or interaction.
	'''

	from pybedtools import BedTool;	
	
	doublebed  = [];
	for c, interval in enumerate(BedTool(path)):
		doublebed.append(interval);
		if(c % 2):
			yield doublebed;
			doublebed = [];
			
			
def generator_seqrecord(paths, ftype):
	'''Generates Bio.SeqRecord.SeqRecord objects from multiple files
	
		paths list: list of paths to files with seqrecords(genbank, fastq, fasta, ets.)
		ftype str: format of files with seqrecords('genbank', 'fastq', 'fasta', ets.)
		
	Yields Bio.SeqRecord.SeqRecord
	'''
	from Bio import SeqIO;	
	for path in paths:
		for seqrecord in SeqIO.parse(path, ftype):
			yield seqrecord;
			
			
def generator_seq(paths, ftype):
	'''Generates sequences from fasta/fastq/genbank files. NOTE: sequences will be converted to upper case
	
		paths list: list of paths to files with seqrecords(genbank, fastq, fasta, ets.)
		ftype str: format of files with seqrecords('genbank', 'fastq', 'fasta', ets.)
		
	Yiels str: nucleotide sequence
	'''
	from Bio import SeqIO;	
	for path in paths:
		for seqrecord in SeqIO.parse(path, ftype):
			yield str(seqrecord.seq.upper());
			
			
			
def generator_mirna(paths, seed_start=1, seed_stop=7):
	'''Generates miRNAs from given fasta files
	
		paths list: list of paths to fasta files with miRNAs
		seed_start int: start postion of mirna seed, 0-based and inclusive
		seed_stop int: stop postion of mirna seed, 0-based and exclusive
		
	Yields nrlbio.Mirna object
	'''
	from Bio import SeqIO;
	from nrlbio.mirna import Mirna
	for path in paths:
		for seqrecord in SeqIO.parse(path, 'fasta'):
			yield Mirna(seqrecord.id, str(seqrecord.seq.upper()), seed_start = seed_start, seed_stop = seed_stop)
			

			
		
	
def grouper(iterable, n):
	it = iter(iterable);
	arr = [];
	while(it):
		for i in range(n):
			try:
				arr.append(next(it));
			except:
				yield arr;
				return
		yield arr;
		arr = []
				
#testing			
if(__name__ == "__main__"):
	#mygen = generator_fastq(sys.argv[1], take = ["id", "seq", "sign", "qual"], shuffle = True)
	#for fq in mygen:
		#print fq;
		#print;
	from conservation import conservation_perfect;	
	mygen = generator_maf(sys.argv[1])
	for fq in mygen:
		print fq;
		print conservation_perfect(fq, "AA")
		print;