#! /usr/bin/python

def shuffle_string(s):
	import random
	l = list(s)
	random.shuffle(l)
	return ''.join(l)

def fasta2dict(path, reverse = False, shuffle = False):
	from Bio import SeqIO;
	"""creates dictionary from fasta file. keys are fasta ids, values are fasta sequences"""
	f=open(path)
	t = SeqIO.to_dict(SeqIO.parse(f,"fasta"))
	f.close()
	
	a = {};
	for k, v in t.iteritems():
		if(shuffle):
			a[k] = shuffle_string(str(v.seq)).upper();
		elif(reverse):
			a[k] = str(v.seq).upper()[::-1];
		else:	
			a[k] = str(v.seq).upper()
	return a;


def counter2string(counter, fraction = False):
	'''produces nice representation of given counter like dictionary'''
	string = ""
	if(fraction):
		t = sum(counter.values())
		for key in sorted(counter.keys()):
			string += "%s\t%1.4f\n" % (key, counter[key]/float(t))		
	else:		
		for key in sorted(counter.keys()):
			string += "%s\t%1.1f\n" % (key, counter[key])
	return string;	
	
def counters2string(*args):
	string = "";
	keys = set();
	for el in args:
		keys.update(el.keys());
	for key in sorted(list(keys)):
		string += str(key) + "\t" + "\t".join([str(x[key]) for x in args]) + "\n";
	return string;	
			
	
def tsv2dict(f, key_collumn, sep = "\t"):
	h = open(f);
	ans = {}
	for line in h:
		arr = line.strip().split(sep);
		ans[arr[key_collumn]] = arr;
	return ans	
	h.close();	
	
def mir2seed(mirdict, start = 1, end = 7):
	from Bio.Seq import reverse_complement;
	sdict = {};
	for k,v in mirdict.items():
		sdict[k] = reverse_complement(v[start: end])
	return sdict;
	
def findsubstring(string, substring, intersect = False):
	positions = [];
	pointer = 0;
	while True:
		pos = string[pointer:].find(substring)
		if(pos == -1):
			return positions;
		else:
			positions.append(pos + pointer);
			if(intersect):
				pointer += pos + 1;
			else:
				pointer += pos + len(substring);
	return positions;
	
	
#print findsubstring("FFFFFFAAAAAFFFFFFFAAAAAFFFFFAAAA", "AAAA", 1)	
	