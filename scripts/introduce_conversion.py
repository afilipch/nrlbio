#! /usr/bin/python
'''script introduces selected conversion into each read of given fastq file. introduces conversion into provided fastq record. For each occurence of nucleotide "from_" at position i generates new fastq where nucleotide at position i is converted to nucleotide "to". For example: convert(fastq("@id1", "TAAAT", "+", "QUALL"), "T", "C") will generates two fastq records: fastq("@id1_T:C:0", "CAAAT", "+", "QUALL"), fastq("@id1_T:C:4", "TAAAC", "+", "QUALL")'''
import argparse;
from nrlbio import generators



parser = argparse.ArgumentParser(description='script introduces selected conversion into each read of given fastq fileintroduces conversion into provided fastq record. For each occurence of nucleotide "from_" at position i generates new fastq where nucleotide at position i is converted to nucleotide "to". For example: convert(fastq("@id1", "TAAAT", "+", "QUALL"), "T", "C") will generates two fastq records: fastq("@id1_T:C:0", "CAAAT", "+", "QUALL"), fastq("@id1_T:C:4", "TAAAC", "+", "QUALL")');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to fastq file");
parser.add_argument('-f', '--from_', nargs = '?', type = str, help = "Conversion from");
parser.add_argument('-t', '--to',  nargs = '?', type = str, help = "Conversion to");
parser.add_argument('-n', '--threads', nargs = '?', default = 4, type = int, help = "number of threads");
args = parser.parse_args();

def convert(fastq, from_, to):
	'''introduces conversion into provided fastq record. For each occurence of nucleotide "from_" at position i generates new fastq where nucleotide at position i is converted to nucleotide "to".
	For example: convert(fastq("@id1", "TAAAT", "+", "QUALL"), "T", "C") will generates two fastq records: fastq("@id1_T:C:0", "CAAAT", "+", "QUALL"), fastq("@id1_T:C:4", "TAAAC", "+", "QUALL")
	
		fastq namedtuple(generators.fastq): initial fastq record, source for converted fastq record
		from_ str: nucleotide to be converted into nucleotide "to"
		to str: nucleotide to be converted in
		
		Yields namedtuple(generators.fastq): fastq record with one conversion, position and type are added to id
	'''	
	for i, n in enumerate(fastq.seq):
		if(n == from_):
			a = list(fastq.seq);
			a[i] = to;
			new_seq = ''.join(a);
			new_id = "%s_%s:%s:%d" % (fastq.id, from_, to, i);
			yield generators.fastq(new_id, new_seq, fastq.sign, fastq.qual)

	


for i, fastq in enumerate(generators.generator_fastq(args.path, take = ["id", "seq", 'sign', 'qual'])):
	for nf in convert(fastq, args.from_, args.to):
		print "\n".join(nf)
