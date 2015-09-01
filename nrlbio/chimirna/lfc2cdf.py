#! /usr/bin/python	
'''produces table of LFC CDF for tables of genes'''
import sys;
import copy;
import os;
import argparse;
from collections import *;
import itertools
from scipy.stats import ks_2samp, mannwhitneyu


parser = argparse.ArgumentParser(description='Script outputs presence of certain binding modes in the interactions');
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
		
		
gene_list = [];	

for path in args.path:
	gene_list.append(get_genes(path));
		
lfc = {};	
f = open(args.lfc);
for l in f:
	a = l.strip().split("\t");
	try:
		lfc[a[0]] = float(a[1]);
	except:
		pass;
f.close();	

lfc_list =[lfc.values()];
for genes in gene_list:
	tlfc = [];
	for g in genes:
		try:
			tlfc.append(lfc[g]);
		except:
			pass;
			#print sys.stderr.write("%s\n" % g);
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
	
for i in range(len(lfc_list)):
	print ('set n%d\tlength %d' % (i, len(lfc_list[i])))
			
#>>> output p-values:
for k1, k2 in itertools.combinations(list(range(len(lfc_list))), 2):
	print "differnce between gene set %d and gene set %d" % (k1+1,k2+1)
	print sys.stderr.write("KS statistics \t%1.3f\tp value\t%.2e\n" % ks_2samp(lfc_list[k1], lfc_list[k2]));
	print


#print sys.stderr.write("KS statistics \t%1.3f\tp value\t%.2e\n" % ks_2samp(lfc_list[1] + lfc_list[2], lfc_list[0] ));



	