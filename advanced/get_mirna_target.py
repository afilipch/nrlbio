#! /usr/lib/python
'''produces old-fashion interaction.bed file of mirna and their targets'''
import argparse;
import os;




#import pysam;
#from nrlbio.filters_for_sam import *
#from nrlbio.chimera import arlist2chimera
#from nrlbio import chimera

parser = argparse.ArgumentParser(description='produces old-fashion interaction.bed file of mirna and their targets');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the new-style interaction file");
parser.add_argument('-m', '--mirna', nargs = '?', required = True, type = str, help = "path to miRNA fasta file");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome ce6|mm9|hg19");
args = parser.parse_args();

exec("from sequence_data.systems import %s as gsys" % args.system);