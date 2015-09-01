import parsing;
from collections import *
from Bio.Seq import reverse_complement

Mapread = namedtuple("Mapped_reads", "chromosome, genome_start, genome_end, read_id, score, strand, read_start, read_end, collapse, conversion, map_quality, cigar, match_part, extrefseq, extupstream, extdownstream")
Chimera = namedtuple("Chimera", "chipart, mapread")


class Cluster(object):
	def __init__(self, mapreads, mid):
		self.mapreads = mapreads;
		self.mid = mid;
		self.start = min([x.genome_start for x in mapreads])
		self.end = max([x.genome_end for x in mapreads])
		self.chromosome = mapreads[0].chromosome
		self.strand = mapreads[0].strand
		self.totalreads = sum([x.collapse for x in mapreads])
		self.indreads = len(mapreads)
		self.indseqs = len(set([x.match_part for x in mapreads]));
		self.score = max([x.score for x in mapreads])
		self.map_quality = max([x.map_quality for x in mapreads])
	#def set_extseq(self, lext, rext, seq):
		#self.extupstream, self.extdownstream, self.extseq = lext, rext, seq;

		
class Interaction(object):
	def __init__(self, cluster, mirid):
		self.cluster = cluster;
		self.mirid = mirid;
		
	def set_mirsseq(self, mirdict, famdict):
		self.mirseq = mirdict.get(self.mirid, "undef");
		if(self.mirseq != "undef"):
			self.seed = reverse_complement(self.mirseq[1:7]);
			return 0;
		else:
			self.seed = reverse_complement(famdict[self.mirid][1:7]);
			self.longseed = famdict[self.mirid]
			return 1;
			
	def extend(self, extu, extd):
		self.extu = extu;
		self.extd = extd;
		self.strand = self.cluster.strand;
		self.chromosome = self.cluster.chromosome;
		
		if(self.strand == "+"):
			self.start = self.cluster.start - extu
			self.end = self.cluster.end + extd
		elif(self.strand == "-"):
			self.start = self.cluster.start - extd
			self.end = self.cluster.end + extu	
		else:
			raise Exception("weird strand!\n")
		
	def set_tseq(self, tseq):
		self.tseq = tseq;
	#def findseed(self):
		#self.spos = parsing.findsubstring(self.cluster.extseq, self.seed);
			
		
def get_interactions(chimeras):
	interactions = [];
	c = 1;
	stripes = defaultdict(list);
	for ch in chimeras:
		stripes[(ch.mapread.chromosome, ch.mapread.strand, ch.chipart.ref_id)].append(ch);
		
	for chrstr, clist in stripes.iteritems():
		clist.sort(key = lambda x: x.mapread.genome_start)
		container = [clist[0]];
		uboundary = clist[0].mapread.genome_end		
		for ch in clist[1:]:
			if(uboundary >= ch.mapread.genome_start):
				container.append(ch);
				uboundary = max(ch.mapread.genome_end, uboundary);
			else:	  
				interactions.append(Interaction(Cluster([x.mapread for x in container], c), container[0].chipart.ref_id));
				container = [ch];
				uboundary = ch.mapread.genome_end
				c += 1;
		interactions.append(Interaction(Cluster([x.mapread for x in container], c), container[0].chipart.ref_id));
		c += 1;
	return interactions;
		
		
def get_clusters(mapreads):
	clusters = [];
	c = 1;
	stripes = defaultdict(list);
	for m in mapreads:
		stripes[(m.chromosome, m.strand)].append(m);

	for chrstr, mlist in stripes.iteritems():
		#print chrstr, len(mlist)
		mlist.sort(key = lambda x: x.genome_start)
		container = [mlist[0]];
		uboundary = mlist[0].genome_end		
		for m in mlist[1:]:
			if(uboundary >= m.genome_start):
				container.append(m);
				uboundary = max(m.genome_end, uboundary);
			else:	  
				clusters.append(Cluster(container, c))
				container = [m];
				uboundary = m.genome_end
				c += 1;
		clusters.append(Cluster(container, c))
		c += 1;
	return clusters;	