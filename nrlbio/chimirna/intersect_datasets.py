#! /usr/bin/python	
'''Script outputs presence of certain binding modes in the interactions'''
import sys;
import os;
import argparse;
from collections import *;


parser = argparse.ArgumentParser(description='Script outputs presence of certain binding modes in the interactions');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (output/interactions.bed) from different dataset");
args = parser.parse_args();





Interaction = namedtuple("Interaction", "chromosome, start, end, iid, score, strand, mirid, mirseq, tseq, map_quality, totalreads, indreads, indseqs, dataset");

class Locus(object):
	def __init__(self, mid, interactions):
		self.mid = mid
		self.interactions = interactions;
		self.start = min([x.start for x in interactions])
		self.end = max([x.end for x in interactions])
		self.score = max([x.score for x in interactions])
		self.chromosome = interactions[0].chromosome;
		self.strand = interactions[0].strand;	
		self.datasets = set([x.dataset for x in self.interactions]);
		self.totalreads = sum([x.totalreads for x in interactions])
		self.indreads = sum([x.indreads for x in interactions])
	
	def represent(self):
		arr = self.chromosome, self.start, self.end, "i%d" % self.mid, self.score, self.strand;
		return "\t".join([str(x) for x in arr])
		
		
		
def get(path, dataset):
	interactions = []
	f = open(path);
	for l in f:
		chromosome, start, end, iid, score, strand, mirid, mirseq, tseq, map_quality, totalreads, indreads, indseqs = l.strip().split("\t")[:13];
		interactions.append(Interaction(chromosome, int(start), int(end), iid, int(score), strand, mirid, mirseq, tseq, int(map_quality), int(totalreads), int(indreads), int(indseqs), dataset));
	f.close()
	return interactions;		
	
def interactions2loci(interactions, overlap = 16):
	loci = [];
	d = defaultdict(list);
	for inter in interactions:
		d[inter.chromosome, inter.strand, inter.mirid].append(inter);
	c=1;	
	for ilist in d.values():
		ilist.sort(key = lambda x: x.start)
		container = [ilist[0]];
		uboundary = ilist[0].end		
		for ch in ilist[1:]:
			if(uboundary >= ch.start + overlap):
				container.append(ch);
				uboundary = max(ch.end, uboundary);
			else:	  
				loci.append(Locus(c, container));
				container = [ch];
				uboundary = ch.end
				c += 1;
		loci.append(Locus(c, container));
		c+=1;
	return loci;
	
	
	
interactions = [];
for path in args.path:
	interactions += get(path, path);
	
loci = interactions2loci(interactions);

intersection = defaultdict(int)
for locus in loci:
	if (locus.indreads > 1):
		intersection[tuple(sorted(locus.datasets))] += 1;
		
for k,v in intersection.iteritems():
	print k, v
			
	
	

