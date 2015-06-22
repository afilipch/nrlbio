#! /usr/lib/python
'''Converts phred scores from one base to another'''
import argparse;



parser = argparse.ArgumentParser(description='Converts integer qualities to phred scores with a given base. or modern sequencing phred base is 33(Illumina 1.8). For previous sequencing versions base is 64 (Starting with Illumina 1.3 and before Illumina 1.8) or 59 (Solexa/Illumina 1.0)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to fastq file");
parser.add_argument('--input_base', nargs = '?', type = str, choices = [33 , 64], help = "phred base of input file");
parser.add_argument('--output_base', nargs = '?', type = str, choices = [33, 64], help = "phred base of output file");
args = parser.parse_args();


def convert_qual(quality, ibase, obase):
	q = [ord(x) - ibase for x in quality] 
	l = [chr(x + obase) for x in q];
	return "".join(l);
	
	
with open(args.path) as f:
	entry = [];
	for l in f:
		entry.append(l.strip())
		if(len(entry)==4):
			print "%s\n%s\n%s" % tuple(entry[:3]);
			print convert_qual(entry[3], args.input_base, args.output_base)
			entry[:] = []