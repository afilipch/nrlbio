#! usr/bin/python
'''script generates stdout report and plots for given fasta file'''
import os;
import sys;
import argparse;

from nrlbio.statistics import fasta
from nrlbio import html




parser = argparse.ArgumentParser(description='script generates stdout report and plots for given fasta/fastq file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file"); 
parser.add_argument('--ftype', nargs = '?', default = 'fasta', choices = ['fasta', 'fastq'], type = str, help = "format of input file. Can be fasta or fastq");
parser.add_argument('--draw', nargs = '?', default = False, const = True, type = bool, help = "if True, scripts also generate plots");
parser.add_argument('--sparse_coefficient', nargs = '?', default = 1, type = int, help = "If sparse_coefficient set to 3, then only each 3rd record will be analyzed, saves time in a case of large files");
args = parser.parse_args();


stat = fasta.Stat.from_file(args.path, args.ftype, sparse_coefficient=args.sparse_coefficient)
print stat


