import re;
from collections import *;


class Sam_record(object):
	def __init__(self, read_id, ref_id, collapse, conversion, start, map_quality, cigar, sequence, score):    
		### set other attributes    
		self.read_id = read_id;
		self.ref_id = ref_id;
		self.collapse = collapse;
		self.conversion = conversion;
		self.ref_start = start; # 0-based
		self.map_quality = map_quality;
		self.cigar = cigar;
		self.sequence = sequence;
		self.score = score;
		self.analize_cigar();
		
		self.refseq = "";
		self.genome_start = 0;
		self.genome_end = 0;
		self.chromosome = "";
		self.strand = "";		
		
	def analize_cigar(self):
		soft_start, soft_end, match, deleted, inserted = (0,0,0,0,0)
		pattern = r"[0-9]+[A-Z]{1}"
		regex = re.compile(pattern);
		clipped = regex.findall(self.cigar);
		if (clipped[0][-1] == "S"):
			soft_start = int(clipped[0][:-1])
		if (clipped[-1][-1] == "S"):
			soft_end = int(clipped[-1][:-1])
		for el in clipped:
			if(el[-1] == "M"):
				match += int(el[:-1])
			elif(el[-1] == "I"):
				inserted += int(el[:-1])  
			elif(el[-1] == "D"):
				deleted += int(el[:-1])      
		self.read_start = soft_start;
		self.read_end = soft_start + match + inserted
		self.ref_end = self.ref_start + match + deleted;
		self.match_part = self.sequence[self.read_start: self.read_end]
		#self.left_part = self.sequence[:soft_start]
		#self.right_part = self.sequence[len(self.sequence) - soft_end:]   
	
	def get_alignment(self):
		alignment = ["",""];
		sp = 0;
		rp = 0;
		pattern = r"[0-9]+[IDM]{1}"
		regex = re.compile(pattern);
		clipped = regex.findall(self.cigar);	
		for el in clipped:
			if(el[-1] == "M"):
				l = int(el[:-1])
				alignment[0] += self.match_part[sp:sp+l]
				alignment[1] += self.refseq[rp:rp+l]
				sp += l
				rp += l
			elif(el[-1] == "I"):
				l = int(el[:-1])  
				alignment[0] += self.match_part[sp:sp+l]
				alignment[1] += "-"*l
				sp += l
			elif(el[-1] == "D"):
				l = int(el[:-1])  
				alignment[0] += "-"*l
				alignment[1] += self.refseq[rp:rp+l]
				rp += l;
		return alignment
		
	def get_conversions(self):
		alignment = self.get_alignment()
		if(self.conversion>-1):
			conv =  [(self.conversion, 'T', "C")];
		else:
			conv = [];
		for i in range(len(alignment[0])):
			if(alignment[0][i] != alignment[1][i]):
				conv.append((i, alignment[1][i], alignment[0][i]));
		return conv;
			
		
def intersect(set1, set2, d = 18):
	'''as input produce the list of tuples: 1st-> sam_record from set1, 2nd -> sam_record from set2'''
	r =[]
	dict1 = defaultdict(list);
	for sr in set1:
		dict1[sr.ref_id].append(sr);
	dict2 = defaultdict(list);
	for sr in set2:
		dict2[sr.ref_id].append(sr);	
		
	for k, lv in dict1.iteritems():
		for v1 in lv:
			for v2 in dict2[k]:
				if(abs(v1.ref_start - v2.ref_start) <= d):
					r.append(v1);
					break;	
	#print len(filter(lambda x: x.score> 34, r))				
	#print len(set([x.read_id for x in filter(lambda x: x.score> 34, r)])) 
	#print Counter([x.score for x in r])
	return r;
	
	
#get_alignment("CACCCCCAAAATTTCGATCTCTCTACTCTTCG", "CACCCCCAAAATTTCGATCTCTCCTACTCTTCG", "22M1D10M")

	
		
		
