#! /usr/bin/python
'''script converts qseq files and sorts them according to TrueSeq barcode present in accompanying files. Ouputs fastq file to a folder provided'''
import sys, os, re, argparse, glob, gzip;
from collections import *;
from Bio import SeqIO;
from Bio.Alphabet import IUPAC
#from Bio.SeqRecord import SeqRecord;
#from Bio.Seq import Seq;



parser = argparse.ArgumentParser(description='script converts qseq files and sorts them according to TrueSeq barcode present in accompanying files');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, default = os.getcwd(), help = "path to the folder with qseq files");
parser.add_argument('-l', '--lane', nargs = '?', required = True, type = int, help = "number of lane");
parser.add_argument('-p', '--pair', nargs = '?', default = 1, type = int, help = "number of pair. Provide in case of paired-end sequencing");
parser.add_argument('-b', '--barcodes', nargs = '?', default = False, type = str, help = "path to the barcodes fasta file. Id of each entry should represent the name of experiment, sequence - barcode sequnce");
parser.add_argument('-o', '--output', nargs = '?', default = "", type = str, help = "path to the folder with results");
parser.add_argument('-u', '--remove_uncalled', nargs = '?', default = False, const = True, type = bool, help = "if set, all reads with uncalled bases will be removed");
parser.add_argument('-gz', '--gzip', nargs = '?', default = False, const = True, type = bool, help = "if set, ouputs gz files");
###optional arguments to be compatible with previous version. One is not encouraged to use them
parser.add_argument('-s', '--barseq', nargs = '+', default = False, type = str, help = "barcode sequence/sequences. Will be ignored if used with \'--barcodes\' argument");
parser.add_argument('-n', '--name', nargs = '+', default = False, type = str, help = "name of output file/files. use with \'--barseq\' argument. Order of file name has to be in corcondance with order of barcodes");
args = parser.parse_args();

class Stat:
	'''class collects statistic of conversion. 
	Cookbook to use. First create an instance, then for each result of process_line() call stat.instance.increment(), then output (str(stat))
	
	Attributes:
	barcodes dict: barcodes used in demultiplexing;
	total int: number of reads processed
	uncalled int: number of reads with uncalled bases removed
	low_quality int: number of reads with low quality removed
	broken_barcodes int: number of reads with undetectable barcode removed
	self.passed dict: Keys are barcode names, Values are numbers of reads found with corresponding barcode 
	'''
	
	def __init__(self, barcodes):
		self.barcodes = barcodes;
		self.total = 0;
		self.uncalled = 0;
		self.low_quality = 0;
		self.broken_barcodes = 0;
		self.passed = defaultdict(int);
		
	def increment(self, a):
		self.total += 1;
		if(len(a) == 2):
			self.passed[a[0]] += 1;
		elif(a == "low_quality"):
			self.low_quality += 1;
		elif(a == "uncalled"):
			self.uncalled += 1;
		else:
			self.broken_barcodes += 1;
			
	def __str__(self):
		s = "total:\t%d\npassed:\t%d\nlow quality:\t%d\nuncalled bases:\t%d\nbroken barcodes:\t%d\n\nbarcodes found\n" % (self.total, sum(self.passed.values()), self.low_quality, self.uncalled, self.broken_barcodes)
		for k, v in self.barcodes.items():
			s += "%s:\t%d\n" % (v, self.passed[v])
		return s;	

def process_line(sequence_line, barcode_line, barcodes, barcode_length, remove_uncalled):
	'''converts qseq lines into fastq entry, amd assignes barcode name to it(demultiplexing). If fails to create fastq entry, reports string explaining the problem
	
		sequence_line str: line from qseq/trueseq sequence file
		barcode_line str: line from qseq/trueseq barcode file
		barcodes dict: Keys:  barcode's sequences, Values: barcode's name which will be used as filenames
		remove_uncalled bool: if True, qseq reads with uncalled bases will bi NOT outputed
	
		Returns tuple or string: 1st element str: barcode's name which will be used as filename, 2nd element str:fastq entry.
		If fails to create fastq entry, reports string explaining the problem
	'''
	sa = sequence_line.strip().split("\t")
	ba = barcode_line.strip().split("\t")
	#check for read quality, 11th field in qseq field is '1' if the read is trustworthy, '0' otherwise
	if(not int(sa[10])):
		return "low_quality"
	#check for uncalled base	
	if(remove_uncalled and "." in sa[8]):
		return "uncalled"
	
	barcode_name = barcodes.get(ba[8][:barcode_length], False)
	#check for presence of sequenced barcode in set of barcodes provided
	if(not barcode_name):
		return ba[8];
		
	#if everything fine we can output the fastq entry with corresponding barcode name
	rid="_".join(sa[0:6])
	seq = sa[8].replace(".", "N")
	qual = sa[9];
	return barcode_name, "@%s\n%s\n+\n%s\n" % (rid, seq, qual);


#parsing barcodes from possible sources into barcode dict. Keys:  barcode's sequences, Values: barcode's name which will be used as filenames
barcodes = OrderedDict();
if(args.barcodes):
	for sq in SeqIO.parse(args.barcodes, "fasta", alphabet = IUPAC.unambiguous_dna):
		barcodes[str(sq.seq)] = sq.id;
elif(args.barseq):
	if(len(args.barseq) == len(args.name)):
		for k, seq in zip(args.name, args.barseq):
			barcodes[seq.upper()] = k;
	else:
		sys.stderr.write('number of items passed to \'--name\' or \'--barseq\' attributes should be equal\n');
		sys.exit();			
else:
	sys.stderr.write('\'--barcodes\' or \'--barseq\' attribute should be provided\n');
	sys.exit();
#infer length of barcodes:
barcode_length = len(barcodes.keys()[0])
	
	

#create statistic object	
stat = Stat(barcodes);	
#compilling all qseq sequence files(trueseq_seq_files) with 	corresponding qseq barcode files (trueseq_barcode_files)
path = os.path.join(os.path.abspath(args.path), '')
trueseq_seq_files = sorted(glob.glob(('%ss_%d_%d_*qseq.txt' % (path, args.lane, args.pair))))
trueseq_barcode_files = sorted(glob.glob(('%ss_%d_%d_*qseq.txt' % (path, args.lane, 2))))
#open all output files:
of = {};
if(args.gzip):
	for barcode_name in barcodes.values():
		of[barcode_name] = gzip.open(os.path.join(args.output, barcode_name + ".fastq.gz"), 'wb')
else:	
	for barcode_name in barcodes.values():
		of[barcode_name] = open(os.path.join(args.output, barcode_name + ".fastq"), 'w')

#run through each line in each file.
for sf, bf in zip(trueseq_seq_files, trueseq_barcode_files):
	sh = open(sf);
	bh = open(bf);
	for sl, bl in zip(sh, bh):
		a = process_line(sl, bl, barcodes, barcode_length, args.remove_uncalled)
		stat.increment(a);
		if(len(a)==2):
			of[a[0]].write(a[1]);
	sh.close();
	bh.close();
#close all output files;
for f in of.values():
	f.close();
print stat;


	

