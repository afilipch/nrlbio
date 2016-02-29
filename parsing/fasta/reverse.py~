#! /usr/bin/python
'''Reverses sequences, provided in fasta format'''
import argparse;

from Bio import SeqIO

parser = argparse.ArgumentParser(description='Reverses sequences, provided in fasta format. NOTE: \'reversed_\' prefix will be added automatically to sequence headers');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to fasta file");
parser.add_argument('--block_size', nargs = '?', type = int, default = 1, help = "If set, blocks of nucleotides are reversed. If set to 2, then dinucleotides will be reversed")
parser.add_argument('-ll', '--line_length', nargs = '?', default = 50, type = int, help = "Max length of line in output fasta file. Makes it more readable");
args = parser.parse_args();


def reverse_fast(seqrecord, line_length):
	line = [];
	for s in seqrecord.seq[::-1]:
		line.append(s);
		if(len(line) == line_length):
			yield "".join(line);
			line[:] = [];
	else:
		yield "".join(line);


def reverse_block(seqrecord, block_size, line_length):
	block = [];
	line = [];
	for s in seqrecord.seq[::-1]:
		block.append(s);
		
		if(len(block)==block_size):
			line.extend(block[::-1])
			block[:] = [];
			
		if(len(line)>=line_length):
			yield "".join(line[:line_length]);
			line = line[line_length:]
			
	else:
		yield "".join(line);
		


for seqrecord in SeqIO.parse(args.path, "fasta"):
	print ">reversed_%s" % seqrecord.id;
	
	if(args.block_size == 1):
		for s in reverse_fast(seqrecord, args.line_length):
			print s;
			
	else:
		for s in reverse_block(seqrecord, args.block_size, args.line_length):
			print s;

		
		