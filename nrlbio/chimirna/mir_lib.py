import re;
import gconst;
import parsing
import os;
from collections import *;
from Bio.Seq import reverse_complement;
from itertools import combinations, product, combinations_with_replacement;

Mismatch = namedtuple("Mismatch", "fr, to, position");
Deletion = namedtuple("Deletion", "fr, position");
Insertion = namedtuple("Insertion", "to, position");# "position" means BEFORE what position in the original sequence insertion happened

class Mismatched(object):
	def __init__(self, initial_sequence, mismatches):
		l = list(initial_sequence);
		#print l
		for m in mismatches:
			l[m.position] = m.to;
		self.seq = "".join(l);	
		self.mismatches = mismatches;
		
class Inserted(object):
	def __init__(self, initial_sequence, insertions):
		l = list(initial_sequence);
		for i in sorted(insertions, key = lambda x: x.position, reverse = True):
			l.insert(i.position, i.to)
		self.seq = "".join(l);	
		self.insertions = insertions;
		
class Deleted(object):
	def __init__(self, initial_sequence, deletions):
		l = list(initial_sequence);
		for d in deletions:
			l[d.position] = ""
		self.seq = "".join(l);	
		self.deletions = deletions;		

class Seed(object):
	def __init__(self, mirseq, mirid, start=1, end = 8):
		self.id = mirid;
		self.mirseq = mirseq;
		self.seed = mirseq[start:end];
		self.match = reverse_complement(self.seed);
		self.mismatched = [];
		self.inserted = [];
		self.deleted = [];
		self.start, self.end = start, end;
		#print self.seed
		#print self.match;
		
	def generate_mismatches(self, upto = 1, start = False, end = False):
		if(start or end):
			match = reverse_complement(self.mirseq[start:end])
		else:
			match = self.match;
			

		level = [];
		for positions in combinations(range(len(match)), upto):
			locker = [match[p] for p in sorted(positions)]
			variants_raw = product("ACTG", repeat = len(positions));
			variants = [];
			for v in variants_raw:
				for i in range(len(v)):
					if(v[i] == locker[i]):
						break;
				else:
					variants.append(v);
			for v in variants:
				mismatches = [];
				for i in range(len(v)):
					mismatches.append(Mismatch(match[positions[i]], v[i], positions[i]))
				level.append(Mismatched(match, mismatches))
		self.mismatched.append(level);
		return None;
		
	def generate_deletions(self, upto = 1):
		for num in range(1,1+upto):
			level = [];
			for positions in combinations(range(1, len(self.match)-1), num):
				deletions = [];
				for p in positions:
					deletions.append(Deletion(self.match[p], p))
				level.append(Deleted(self.match, deletions))
			#for d in level:
				#print d.seq, d.deletions;
			#print "\n"	
			self.deleted.append(level);
		return None;	
		
	def generate_insertions(self, upto = 1):
		for num in range(1,1+upto):
			level = [];
			for positions in combinations_with_replacement(range(1, len(self.match)), num):
				variants = product("ACTG", repeat = len(positions));
				for v in variants:
					insertions = [];
					for i in range(len(v)):
						insertions.append(Insertion(v[i], positions[i]))
					level.append(Inserted(self.match, insertions))
			#for d in level:
				#print d.seq, d.insertions;
			#print "\n"	
			self.inserted.append(level);
		return None;
		
	def generate_default_variants(self):
		self.generate_mismatches(upto = 1);
		self.generate_mismatches(upto = 2, start = 1, end = 8);
		self.generate_deletions(upto = 1);
		self.generate_insertions(upto = 1);
		
	def associated_values(self, mtype):
		def get_ass_list(mtype, reflist):
			r = [];
			for i, el1 in enumerate(reflist):
				r.append([]);
				for el2 in el1:
					r[i].append(mtype())
			return r		
			
		self.ass_match = mtype();
		self.ass_mismatched = get_ass_list(mtype, self.mismatched);
		self.ass_inserted = get_ass_list(mtype, self.inserted);
		self.ass_deleted = get_ass_list(mtype, self.deleted);
		
	def find_once(self, seq):
		'''hierarchically searches for match or its variations inside the sequence provided'''
		ans = [];
		
		if(self.match in seq):
			ans.append(self.match);
			return "match", ans;
			
		for mismatched in self.mismatched[0]:
			if(mismatched.seq in seq):
				ans.append(mismatched);
		if(ans):
			return "mismatch", ans
			
		for inserted in self.inserted[0]:
			if(inserted.seq in seq):
				ans.append(inserted);
		if(ans):
			return "insertion", ans		
			
		if(len(self.mismatched) > 1):	
			for mismatched in self.mismatched[1]:
				if(mismatched.seq in seq):
					ans.append(mismatched);
			if(ans):
				return "mismatch2", ans				
			
		for deleted in self.deleted[0]:
			if(deleted.seq in seq):
				ans.append(deleted);
		if(ans):
			return "deletion", ans		
			
				
		return "no_match", ans		
				
	def find_everything(self, seq):
		'''finds every possible presence of match or it variations inside the sequence provided'''
		ans = {'match': [], 'mismatch': [], 'insertion': [], 'deletion': []};

		if(self.match in seq):
			ans['match'].append(self.match);
			
		for mismatched in self.mismatched[0]:
			if(mismatched.seq in seq):
				ans['mismatch'].append(mismatched);

		for inserted in self.inserted[0]:
			if(inserted.seq in seq):
				ans['insertion'].append(inserted);	
			
		for deleted in self.deleted[0]:
			if(deleted.seq in seq):
				ans['deletion'].append(deleted);
				
		return ans;		
		
def get_seeds (mirdict, start=1, end=8):
	seeds = {};
	for mirid, mirseq in mirdict.items():
		seeds[mirid] = Seed(mirseq, mirid, start, end);
		seeds[mirid].generate_default_variants();
		#seeds[mirid].generate_deletions(upto = 1);
		#seeds[mirid].generate_insertions(upto = 1);
	return seeds;	
			
			
		
if __name__ == "__main__":					
	seed = Seed("TGAGGTAGTAGGTTGTATAGTT", "cel-let-7")			
	#seed.generate_mismatches(2);
	#seed.generate_deletions(2);
	seed.generate_insertions(2);
	print seed.match
		
		


def mir2seed(mirdict, start = 1, end = 7):
	sdict = {};
	for k,v in mirdict.items():
		sdict[k] = reverse_complement(v[start: end])
	return sdict;
	
def seed2mirs(mirdict, start = 1, end = 7):	
	sdict = defaultdict(list);
	for k,v in mirdict.items():
		sdict[reverse_complement(v[start: end])].append(k);
	return sdict;	

	
	
def mm28(seq, mir):
	#if(len(mir)> 7):
	s27 = reverse_complement(mir[1:7])
	s38 = reverse_complement(mir[2:8])
	a = reverse_complement(mir[1:8])
	#elif(len(mir) == 7):
		#s27 = reverse_complement(mir[:6])
		#s38 = reverse_complement(mir[1:7])
		#a = reverse_complement(mir);
		
	mm37 = set();
	for i in range(1,6):
		for n in "ACTG":
			mm37.add(a[:i] + n + a[i+1:]);
	order = [s27, s38] + list(mm37);
	for el in order:
		if(el in seq):
			return seq.rfind(el);
	return -1;

	
def modes(seq, mir):
	#if(len(mir)> 7):
	s27 = reverse_complement(mir[1:7])
	s38 = reverse_complement(mir[2:8])
	a = reverse_complement(mir[1:8])
	#elif(len(mir) == 7):
		#s27 = reverse_complement(mir[:6])
		#s38 = reverse_complement(mir[1:7])
		#a = reverse_complement(mir);	
				
	mm28 = set();
	for i in range(1,6):
		for n in "ACTG":
			mm28.add(a[:i] + n + a[i+1:]);
	
	wobble = set()
	for i in range(1,6):
		if(a[i] == "C"):
			wobble.add(a[:i] + "T" + a[i+1:]);
		elif(a[i] == "A"):
			wobble.add(a[:i] + "G" + a[i+1:]);
			
	ins28 = set();# target bulge
	for i in range(1,6):
		for n in "ACTG":
			ins28.add(a[:i] + n + a[i:]);				
			
	del28 = set();# miRNA bulge
	for i in range(1,6):
		del28.add(a[:i] + a[i+1:]);			
			
			
	order = [s27], [s38] , wobble, mm28, ins28, del28;
	for i in range(len(order)):
		for s in order[i]:
			if(s in seq):
				return i;
	return -1;		
	

	
	
	
pattern1 = r"\w+-(\w+)-(\d+.*)"
regex1 = re.compile(pattern1);  

def longname(mirids):
	result = "";
	mirs = set();
	others = set();
  
	for m in mirids:  
		temp = regex.search(m).groups()
		if(temp[0] == "miR" or temp[0] == "mir"):
			mirs.add(temp[1]);
		else:
			others.add("-".join(temp))
	if(others):
		result += ",".join(others) + "\\n"
	elif(mirs):
		result += "miR(" + ",".join(mirs) + ")\\n"
	return result 
  
  
pattern2 = r"\w+-(\w+)-(\d+)"
regex2 = re.compile(pattern2);  
  
def shortname(mirids):
	result = "";
	mirs = set();
	others = set();
	for m in mirids: 
		temp = regex2.search(m).groups()
		if(temp[0] == "miR" or temp[0] == "mir"):
			mirs.add(temp[1]);
		else:
			others.add("-".join(temp))
	if(others):
		return sorted(others)[0] + "\\nfamily"
	else:
		return "miR-" + str(sorted([int(x) for x in mirs])[0]) + "\\nfamily"
		
		
def seeds2name(system):
	mirdict = parsing.fasta2dict(gconst.system2mir[system]);
	s2m = seed2mirs(mirdict);
	s2n = {};
	for k, v in s2m.iteritems():
		s2n[k] = shortname(v);
	return s2n;	
	
def read_families(path):
	'''creates fam2mirs dictionary (keys: fam ids, values sets of miRNAs' ids)
	path String: path to fam_ids.tsv (left/fam_ids.tsv) file;
	'''
	fm = defaultdict(set);
	f = open(path);
	for l in f:
		a = l.strip().split("\t");
		fm[a[0]].update(a[1].split(","));
	f.close();
	return fm;	
	
def get_common_part(seqs, fr = 0, tolerate_first_mismatch = -1):
	'''gets the longest common onset (5' part of mirna) for all given sequences'''
	limit =  min([len(x) for x in seqs])
	ambiguos = 0;
	for i in range(fr, limit):
		if(len(set([x[i] for x in seqs])) > 1): 
			if ((i > tolerate_first_mismatch)):	
				return 'X'*(fr+ambiguos) + seqs[0][fr+ambiguos:i], i;
			else:
				ambiguos = i; 
	return 	'X'*(fr+ambiguos) + seqs[0][fr+ambiguos:limit], limit;
	