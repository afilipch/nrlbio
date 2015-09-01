#! /usr/bin/python	
'''library supporting processing chiparts. Chipart is a string in form: read_id, ref_id, read_start, read_end, read_length, ref_start, ref_end, ref_length, cut_left, cut_right, conversions, weight, score. Chipart represent found part of reference sequence(miRNA) in read sequence'''
from collections import *;

Chipart = namedtuple("Chipart", "read_id, ref_id, read_start, read_end, read_length, ref_start, ref_end, ref_length, cut_left, cut_right, conversions, weight, score");

def read(paths):
	chiparts = [];
	for path in paths:
		f = open(path);
		for line in f:
			line = line.strip();
			read_id, ref_id, read_start, read_end, read_length, ref_start, ref_end, ref_length, cut_left, cut_right, conversions, weight, score = line.split("\t");
			chiparts.append(Chipart(read_id, ref_id, int(read_start), int(read_end), int(read_length), int(ref_start), int(ref_end), int(ref_length), cut_left, cut_right, conversions.split("|"), float(weight), float(score)));
		f.close();
	return chiparts;
	
def repr_chipart(chipart):
	return "\t".join([str(eval("chipart." + key)) for key in ["read_id", "ref_id", "read_start", "read_end", "read_length", "ref_start", "ref_end", "ref_length", "cut_left", "cut_right", "conversions", "weight", "score"]]);
		
	
def get_weighted(chiparts, key):
	ans = defaultdict(float);
	for c in chiparts:
		value = eval("c." + key);
		ans[value] += c.weight;
	return ans;	
	
	
def get_stats(chiparts):
	stat = {};
	for key in ["read_start", "read_end", "ref_start", "ref_end", "cut_right", "cut_left"]:
		stat[key] = get_weighted(chiparts, key)
		
	stat["cut_length"] = defaultdict(float);
	stat["right_length"] = defaultdict(float);
	for c in chiparts:
		stat["cut_length"][c.ref_length - c.ref_end] += c.weight;
		stat["right_length"][c.read_length - c.read_end] += c.weight;
		
	stat["cut_right_after"] = get_weighted(chiparts, 'cut_right[0]') 
	stat["cut_left_after"] = get_weighted(chiparts, 'cut_left[0]') 
	
	stat["num_conversion"] = defaultdict(float);
	stat["pos_conversion"] = defaultdict(float);
	stat["type_conversion"] = defaultdict(float);
	for c in chiparts:
		conv = c.conversions
		if(conv[0] == "-1,N,N"):
			stat["num_conversion"][0] += c.weight;
		else:
			stat["num_conversion"][len(conv)] += c.weight;
		for a in conv:
			b = a.split(",")
			stat["pos_conversion"][int(b[0])] += c.weight;
			stat["type_conversion"]["".join(b[1:3])] += c.weight;						
	return stat;
	
def get_starts(chiparts):
	ans = defaultdict(float);
	for c in chiparts:
		ans[(c.read_start, c.ref_start)] += c.weight;
	return ans;	
	
def repr_stats(stat1, stat2, path, attribute):
	f = open(path, 'w');
	s1, s2 = stat1[attribute], stat2[attribute]
	w1 = sum(s1.values());
	w2 = sum(s2.values())+0.0001;
	for k in sorted(set(s1.keys() + s2.keys())):
		f.write("%s\t%1.3f\t%1.3f\n" % (str(k), s1[k]/w1, s2[k]/w2))
	f.close();	
	
def repr_starts(s1, s2, path):
	f = open(path, 'w');
	w1 = sum(s1.values());
	w2 = sum(s2.values())+0.0001;
	for k in sorted(set(s1.keys() + s2.keys()), key = lambda x: x[0] + x[1]):
		f.write("%d in read %d in ref\t%1.3f\t%1.3f\n" % (k[0], k[1], s1[k]/w1, s2[k]/w2))
	f.close();		

	
	
	