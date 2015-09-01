#! /usr/bin/python	
'''compare pairwise log fold change of transcripts targeted by different miRNA'''
import sys;
import copy;
import os;
import argparse;
from collections import *;
from scipy.stats import ks_2samp, mannwhitneyu
import itertools


parser = argparse.ArgumentParser(description='compare pairwise log fold change of transcripts targeted by different miRNA');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to gene names tables");
parser.add_argument('--lfc', nargs = '?', type = str, help = "path to the table which connects gene names to LFC")
parser.add_argument('-o', '--output', nargs = '?', type = str, help = "name of the output")
args = parser.parse_args();

def get_genes(path):
	genes = set();
	f = open(path);
	for l in f:
		genes.add(l.strip());
	f.close()	
	return genes;	
	
#>>>>>pairwise comparison
def compare(genes1, genes2, lfc):
	i1 = genes1 - genes2;
	i2 = genes2 - genes1;
	lfc1, lfc2 = [],[];
	for g in i1:
		try:
			lfc1.append(lfc[g]);
		except:
			pass;		
	for g in i2:
		try:
			lfc2.append(lfc[g]);
		except:
			pass;
	print len(lfc1)
	print len(lfc2)
	a = [sum(lfc1)/float(len(lfc1)), sum(lfc2)/float(len(lfc2))] + list(mannwhitneyu(lfc1, lfc2))		
	return "mean1 %1.4f\tmean2 %1.4f\nU Mann-Whitney statistics%1.4f\tp value\t%1.9f\n" % (a[0],a[1],a[2],a[3])		
	
#>>> get lfc dict:	key -> gene name, value -> log fold change
lfc = {};	
f = open(args.lfc);
for l in f:
	a = l.strip().split("\t");
	try:
		lfc[a[0]] = float(a[1]);
	except:
		pass;
f.close();	
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>	
		
#>>> get dictionary: key -> name of gene set, value -> gene set		
gene_dict = OrderedDict();
for path in args.path:
	gene_dict[path] = get_genes(path);
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>	


	
	
for k1, k2 in itertools.combinations(gene_dict, 2):
	print "differnce between gene set %s and gene set %s" % (k1,k2)
	print compare(gene_dict[k1], gene_dict[k2], lfc)
	print
	
	



lfc_list =[lfc.values()];
for p, genes in gene_dict.iteritems():
	tlfc = [];
	for g in genes:
		try:
			tlfc.append(lfc[g]);
		except:
			print sys.stderr.write("%s\n" % g);
	lfc_list.append(copy.copy(tlfc));
	
##output
o = open(args.output, 'w')
for i in range(len(lfc)):
	a = []
	for el in lfc_list:
		try:
			v = str(el[i])
		except:
			v = " "			
		a.append(v);
	o.write("\t".join(a) + "\n")
o.close()	
			
			
		






	