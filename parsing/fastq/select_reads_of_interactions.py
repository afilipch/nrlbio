#! /usr/bin/python
'''Selects those reads from the fastq file which form interactions/circles/clusters of interest''' 
import sys;
import argparse




parser = argparse.ArgumentParser(description='Selects those reads from the fastq file which form interactions/circles/clusters of interest');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the sequences, fastq format");
parser.add_argument('--table', nargs = '?', required = True, type = str, help = "Path to chiflex read2iid file, which connects read ids to interactions")
args = parser.parse_args();



selection = set()
with open(args.table) as f:
	for l in f:
		selection.add(l.strip().split("\t")[1]);
		

counter = 0;
with open(args.path) as f:
	entry = [];
	for l in f:
		if(len(entry)==4):
			name = entry[0][1:].split(" ")[0]
			if(name in selection):
				print "\n".join(entry);
			else:
				pass;
			entry[:] = [];
			counter+=1
			if(not counter % 1000000):
				sys.stderr.write("reads processed: %d\n" % counter);
			
		entry.append(l.strip())	