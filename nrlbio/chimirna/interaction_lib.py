#! /usr/bin/python
from collections import *
import random;
import miscellaneous;
from Bio.Seq import reverse_complement;


Interaction = namedtuple("Interaction", "chromosome, start, end, iid, score, strand, mirid, mirseq, tseq, map_quality, totalreads, indreads, indseqs");

		

def get(path, undef = False):
	interactions = []
	f = open(path);
	for l in f:
		chromosome, start, end, iid, score, strand, mirid, mirseq, tseq, map_quality, totalreads, indreads, indseqs = l.strip().split("\t")[:13];
		if(undef or len(mirid.split(",")) == 1):
			interactions.append(Interaction(chromosome, int(start), int(end), iid, int(score), strand, mirid, mirseq, tseq, int(map_quality), int(totalreads), int(indreads), int(indseqs)));
	f.close()
	return interactions;
	
def shufseq(interactions):
	r = [];
	for inter in interactions:
		chromosome, start, end, iid, score, strand, mirid, mirseq, tseq, map_quality, totalreads, indreads, indseqs = list(inter);
		new_tseq = miscellaneous.shuffle_str(tseq);
		r.append(Interaction(chromosome, start, end, iid, score, strand, mirid, mirseq, new_tseq, map_quality, totalreads, indreads, indseqs))
	return r;
	
def shufpairs(interactions):
	r = [];
	ms = [(x.mirid, x.mirseq) for x in interactions];
	for inter in interactions:
		index = random.randint(0,len(ms) - 1);
		m = ms.pop(index);
		chromosome, start, end, iid, score, strand, mirid, mirseq, tseq, map_quality, totalreads, indreads, indseqs = list(inter);
		r.append(Interaction(chromosome, start, end, iid, score, strand, m[0], m[1], tseq, map_quality, totalreads, indreads, indseqs))
	return r;
	
def int2reads(path):
	'''creates int2reads dictionary (keys: interaction ids, values sets of reads' ids)
	path String: path to read2int.tsv file;
	'''
	ir = defaultdict(set);
	f = open(path);
	for l in f:
		a = l.strip().split("\t");
		ir[a[1]].add(a[0]);
	f.close();
	return ir;
	
	
	
class Locus(object):
	def __init__(self, mid, interactions):
		self.mid = mid
		self.interactions = interactions;
		self.start = min([x.start for x in interactions])
		self.end = max([x.end for x in interactions])
		self.score = max([x.score for x in interactions])
		self.chromosome = interactions[0].chromosome;
		self.strand = interactions[0].strand;
		
	def multitargeting(self):	
		self.mirids = set();
		for inter in self.interactions:
			self.mirids.update(inter.mirid.split(","));
		self.seeds = set([reverse_complement(x.mirseq[1:7]) for x in self.interactions])	
		self.s2m = defaultdict(set);
		for inter in self.interactions:
			self.s2m[reverse_complement(inter.mirseq[1:7])].add(tuple(inter.mirid));
		self.onehit = len(self.mirids) == 1 and self.interactions[0].indreads == 1;
		self.onemir = len(self.mirids) == 1;
		self.onefam = False
		for k,v in self.s2m.iteritems():
			if(len(v) > 1 and not bool(set(list(v)[0]).intersection(v))):
				self.onefam = True;
				break;
		self.diffam = len(self.seeds) > 1		
	def represent(self):
		arr = self.chromosome, self.start, self.end, "i%d" % self.mid, self.score, self.strand, "|".join([x.mirid for x in self.interactions]);
		return "\t".join([str(x) for x in arr])
		
		
		
		
	
def interactions2loci(interactions, overlap = 18):
	loci = [];
	d = defaultdict(list);
	for inter in interactions:
		d[inter.chromosome, inter.strand].append(inter);
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
	
	
def merge(interactions, gsys, int2read, name, overlap = 18):
	
	int2newid = defaultdict(set);
	
	def container2interaction(interactions, mid):
		iid = mid
		start = min([x.start for x in interactions])
		end = max([x.end for x in interactions])
		chromosome = interactions[0].chromosome;
		strand = interactions[0].strand;	
		score = max([x.score for x in interactions])
		mirid = interactions[0].mirid;
		mirseq = interactions[0].mirseq;
		tseq = gsys.genome.get_oriented(chromosome, start, end, strand).upper()
		
		map_quality = max([x.map_quality for x in interactions])
		totalreads = sum([x.totalreads for x in interactions])
		indreads = sum([x.indreads for x in interactions])
		indseqs	= 	max([x.indseqs for x in interactions])
		
		return Interaction(chromosome, start, end, iid, score, strand, mirid, mirseq, tseq, map_quality, totalreads, indreads, indseqs);	
		
	
	
	r = [];
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
				r.append(container2interaction(container, '%s%05d' % (name, c)));
				for inter in container:
					int2newid['%s%05d' % (name, c)].update(int2read[inter.iid])
				container = [ch];
				uboundary = ch.end
				c += 1;
		r.append(container2interaction(container, '%s%05d' % (name, c)));
		for inter in container:
			int2newid['%s%05d' % (name, c)].update(int2read[inter.iid])		
		c+=1
	return r, int2newid;	
	
	
	

	
	
	