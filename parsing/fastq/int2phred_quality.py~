#! /usr/lib/python
'''Converts integer qualities to phred scores with a given base'''
import argparse;



parser = argparse.ArgumentParser(description='Converts integer qualities to phred scores with a given base');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to adapters");
parser.add_argument('--phred_base',  default = 33,  nargs = '?', type = int, help = "Phred base. For modern sequencing phred base is 33(Illumina 1.8). For previous sequencing versions base is 64 (Starting with Illumina 1.3 and before Illumina 1.8) or 59 (Solexa/Illumina 1.0)");
args = parser.parse_args();


def convert_qual(quality, base):
	l = [chr(int(x) + base + 5) for x in quality.split()];
	return "".join(l);
	
	
with open(args.path) as f:
	entry = [];
	for l in f:
		entry.append(l.strip())
		if(len(entry)==4):
			print "%s\n%s\n%s" % tuple(entry[:3]);
			print convert_qual(entry[3], args.phred_base)
			entry[:] = [];


