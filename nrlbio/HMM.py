# /usr/bin/python
'''Library of classes and functions related to Markov Models'''

import sys
import random
import itertools
from collections import defaultdict, deque
import copy

import numpy as np

from nrlbio.random_extension import weighted_choice_fast, weighted2interval
from nrlbio.statistics.functions import wilson_score_interval_binomial
from nrlbio import numerictools

conf = {"slide_size": 1000, 'modeldiff': 0.01, 'maxlook': 100, 'jump': 300, 'jumpdiff': 0.009, 'jumplook': 400}

def find_switch(string, growing, sliding, curpos, slide_size=1200, modeldiff=0.01, maxlook=100, jump=300, jumpdiff=0.009, jumplook=400, **kwargs):
	'''Function tries to find an exact position of a swith between different markov models on a given string. For example, for "CCCCCCCCCCCCCCCCCCCCCCCC|TTTTTTTTTTTTTTTTTTTTTTTT" two different markov models should be generated with a switch exactly at "|". NOTE: Solution is merely suboptimal, since additional heuristics is used to speed up an algorithm. Output also depends on heuristic's parameters, so they shoud be either kept default or changed consciously.
	
		Input arguments:	
	
		string str: string to generate separate markov models on
		growing GrowingMarkovChain: markov model, which tries iteratively extend along the string. The switch between it and sliding model has to be found. Shoud be initialized prior call of this function, on a string left to the "string" provided as input.
		sliding SlidingMarkovChain: markov model, which slides "slide_size" chars ahead the growing one in order to calculate a difference between them stepwise. The switch between it and growing model has to be found. Shoud be initialized prior call of this function, on first "slide_size" chars of "string"
		curpos int: current position of a "string" on a global sequence. Global sequence iteratively switches between models, each switch corresponds to a single call of this function with a "string" being substring of a global sequence.
		
		Parameters:
		
		slide_size int: length of sliding model. The bigger the length the faster the algorithm but less precise.
		modeldiff float: minimal difference between growing and sliding windows to stop search and output switch position and growing model. The bigger the modeldiff, the fewer models(but more different from each other) are generated.
		maxlook int: we don't wnat to find a switch at a first position with difference between models more than modeldiff. To prevent it we have to look forward to find the true global optimum(position with the largest difference). The bigger it is, the more precise the algorithm but also more slow.
		jump int: If the difference between growing and sliding model is reasonbly low, it makes sense not to iterate symbol by symbol, but jump "jump" symbols further. The bigger the jump, the faster the algorithm but less precise.
		jumpdiff float: If the difference between growing and sliding model is reasonbly low, it makes sense not to iterate symbol by symbol, but jump "jump" symbols further if difference is less than jumpdiff. The bigger the jumpdiff, the faster the algorithm but less precise.
		jumplook int: jump also can happend if there is no jump or model switch happens for [jumplook] number of iterations. The bigger the jumplook, the faster the algorithm but less precise
		
	Returns tuple:
		growing GrowingMarkovChain: completed growing model, may be stored somewhere for furhter use
		ngrowing GrowingMarkovChain: growing model required for a following iteration of the function
		nsliding SlidingMarkovChain: sliding model required for a following iteration of the function
		position int: position of a model switch on a string
	'''
	
	maxscore = 0;
	lookforward = 0;
		
	
	for i, (ls, rs) in enumerate(zip(string, string[slide_size:])):
		diff = markov_difference(growing, sliding);
		
		if(lookforward==maxlook and maxscore>modelfiff):
			growing = growing.shrink(string[i-lookforward-growing.order*2:i])
			ngrowing = GrowingMarkovChain.from_string(string[i-lookforward: i+slide_size-lookforward], args.order)
			nsliding = SlidingMarkovChain.from_string(string[i+slide_size-lookforward: i+slide_size*2-lookforward], args.order)
			
			sys.stderr.write("switch to a new model happens at position (%d), cause the difference (%1.5f) is more than max difference %1.5f and lookforward is more or equal than %d\n" % 
			(i+curpos-lookforward, maxscore, modelfiff, lookforward))
			return growing, ngrowing, nsliding, i-lookforward+slide_size;
			
		elif(diff<jumpdiff or lookforward>maxlook*5):
			growing.grow_long(string[i: i+jump])
			sliding = SlidingMarkovChain.from_string(string[i+jump: i+slide_size+jump], sliding.order)
			sys.stderr.write("jump of length (%d) happens at position (%d), cause the difference (%1.5f) is less than jump difference %1.5f\n" % (jump, i+curpos, diff, jumpdiff))
			return None, growing, sliding, i+jump;
			
		else:
			growing.grow(ls);
			sliding.slide(rs);	
			if(diff>maxscore):
				lookforward = 0;
				maxscore = diff;
			else:
				lookforward +=1;
	else:
		growing.add(sliding)
		return growing, None, None, 0




class MarkovException(Exception):
	pass;	
	
	

class MarkovChain(object):
	'''Class provides functionality of the basic markov chain. Next state depends on previous one only (transition matrix is 2-dimensional)
	
	Attributes:
		transitions dict of dicts of integers: probabilities(values) of the switch from the previous state (1st layer key) to the next one(2nd layer first list)
		initiation dict of integers: probabilities(values) to start sequence with particular state(key)
		order int: defines the length of the state(for example markov chain with order of two will consider states of the sequence "AGAT" as "AG", "GA", "AT")
		length int: if markov model was generated from sequence then length of the sequence is stored in this attribute
	'''	
	def __init__(self, transitions, initiation, order, seq=None):
		self.transitions = defaultdict(lambda: defaultdict(int), transitions);
		self.initiation = defaultdict(int, initiation);
		self.order = order
		if(seq):
			self.length = len(seq)
		
		self.states = set(initiation.keys() + transitions.keys());
		for d in transitions.values():
			self.states.update(d.keys());
		self._fix_integrity()
		
		self.emsupport = defaultdict(int, dict([(x[0], sum(x[1].values())) for x in self.transitions.items()]))
		self.support = sum(self.emsupport.values())
		
		
	def _fix_integrity(self):
		'''All possible states should get transition and emission pseudocounts'''
		for state in self.states:
			if(state not in self.transitions.keys()):
				self.transitions[state] = {};
				
		for k, d in	self.transitions.iteritems():
			for state in self.states:
				if(state not in d.keys()):
					self.transitions[k][state]=1;
				
			
			
	@classmethod
	def from_string(cls, seq, order):
		'''Creates MarkovChain instance on basis of given sequence
		
			sequence str: sequence to build model on
			order int: defines the length of the state(for example markov chain with order of two will consider states of the sequence "AGAT" as "AG", "GA", "AT")
		'''
		qlen = order*2
		if(len(seq)<=qlen):
			raise MarkovException("it is impossible to generate Markov Model of order %d from the sequence of length(%d) less then double order\n" % (order, len(seq)))
		
		else:
			transitions = defaultdict(lambda: defaultdict(int));
			initiation = defaultdict(int)
			ss = seq[:qlen]
			
			for s in seq[qlen:]:
				transitions[ss[:order]][ss[order:]]+=1
				ss = "".join((ss[1:], s))
			else:
				transitions[ss[:order]][ss[order:]]+=1
				
			for k, d in transitions.iteritems():
				initiation[k] = sum(d.values())
			
		return cls(dict(transitions), dict(initiation), order, seq);
		
		
	def add(self, other):
		for state, d in other.transitions.items():
			for k, v in d.items():
				self.transitions[state][k] += v;
				
		for state, v in other.initiation.items():
			self.initiation[state] += v;
			
		if(self.length and other.length):
			self.length += other.length;
			
		self.states.update(other.states);
		
		self._fix_integrity()
		
		for state, v in other.emsupport.items():
			self.emsupport[state] += v;
		
		self.support += other.support
		
		
	def generate_string(self, length=None, chunk_size=None):
		'''Generates string on basis of the initiation and transition's probabilities'''
		if(length is None):
			length = self.length 
			
		tr = {};	
		for k, d in self.transitions.iteritems():
			tr[k] = weighted2interval(d.items());			
		init = weighted2interval(self.initiation.items())			
				
		last = weighted_choice_fast(*init);
		l = [last]
		curlength = self.order
		
		for _ in range(length/self.order):
			last = weighted_choice_fast(*tr[last]);
			l.append(last);
			curlength+=self.order
			if(curlength>=chunk_size):
				s = "".join(l);
				overhang = s[chunk_size:]
				l = [overhang]
				curlength = len(overhang);
				yield s[:chunk_size];
		else:
			s = "".join(l);
			nl = (length/self.order+1)*self.order - length
			if((len(s)-nl)>0):
				yield s[:len(s)-nl]
				
				
				
class SlidingMarkovChain(MarkovChain):
	'''Class provides functionality of the basic markov chain. Next state depends on previous one only (transition matrix is 2-dimensional)
	
	Attributes:
		transitions dict of dicts of integers: probabilities(values) of the switch from the previous state (1st layer key) to the next one(2nd layer first list)
		initiation dict of integers: probabilities(values) to start sequence with particular state(key)
		order int: defines the length of the state(for example markov chain with order of two will consider states of the sequence "AGAT" as "AG", "GA", "AT")
		seq ...
	'''		
	def __init__(self, transitions, initiation, order, seq):
		super(SlidingMarkovChain, self).__init__(transitions, initiation, order, seq)
		self.seq = seq;
		
		
	def slide(self, rsymbol):
		k = self.order;
		nseq = "".join((self.seq[1:], rsymbol));
		
		self.emsupport[self.seq[:k]] -= 1
		self.emsupport[nseq[-k:]] += 1
		
		self.transitions[self.seq[:k]][self.seq[k:k*2]] -= 1
		self.transitions[nseq[-k*2:-k]][nseq[-k:]] += 1
		
		self.seq = nseq;
		
	#def slide_long(self, seq):
		#l = length(seq)
		#k = self.order
		#nseq = nseq = "".join((self.seq[l:], seq));
		
		#for s in seq[:l+k]:
			

class GrowingMarkovChain(MarkovChain):
	'''Class provides functionality of the basic markov chain. Next state depends on previous one only (transition matrix is 2-dimensional)
	
	Attributes:
		transitions dict of dicts of integers: probabilities(values) of the switch from the previous state (1st layer key) to the next one(2nd layer first list)
		initiation dict of integers: probabilities(values) to start sequence with particular state(key)
		order int: defines the length of the state(for example markov chain with order of two will consider states of the sequence "AGAT" as "AG", "GA", "AT")
		seq ...
	'''		
	def __init__(self, transitions, initiation, order, seq):
		super(GrowingMarkovChain, self).__init__(transitions, initiation, order, seq)
		self.last = seq[-order*2:];
		
	def grow(self, rsymbol):
		k = self.order;
		nlast = "".join((self.last[1:], rsymbol));
		
		self.emsupport[nlast[-k:]] += 1
		self.transitions[nlast[-k*2:-k]][nlast[-k:]] += 1
		self.support+=1
		
		self.nlast = nlast;
		
	def grow_long(self, seq):
		for rsymbol in seq:
			self.grow(rsymbol);
			
			
	def shrink(self, seq):
		k = self.order;
		self.nlast = seq[:k*2];
		nlast = self.nlast
		
		for s in seq[k*2:]:
			nlast = "".join((nlast[1:], s));
			self.emsupport[nlast[-k:]] -= 1
			self.transitions[nlast[-k*2:-k]][nlast[-k:]] -= 1
			self.support-=1
			
			
		
			
			
		
		
class MultiMarkov(object):
	'''Set of connected(through stochastic switches) Markov Models. May be useful for sequences with regions of vastly different transitions' probabilities.
	NOTE: all the models correspond to ONE sequence
	
	Attributes:
		models list of markov models: list of Markov Models, each representing specific region of the sequence. NOTE: all the models correspond to ONE sequence
		switches list of floats: list of probabilities(in interval format) for selection of the models
	'''
	
	
	def __init__(self, models, switches=None):
		self.models = models;
		if(switches):
			self.switches;
		else:
			self.switches = np.linspace(0,1,num=len(models)+1)[1:]

		
	@classmethod	
	def from_string(cls, seq, order, window_size, maxdiff=0.0001):
		models = [];
		l = [];
		for s in seq:
			l.append(s);
			if(len(l)>=window_size):
				mm = MarkovChain.from_string("".join(l), order);
				if(models):
					closest, min_distance = models[0], markov_difference(mm, models[0])
					print "initial distance: % 1.5f\n" % min_distance;
					for cmm in models[1:]:
						distance = markov_difference(mm, cmm);
						print "distance: % 1.5f" % distance;
						if(distance<min_distance):
							min_distance=distance
							closest=cmm;
					print "minimal distance: % 1.5f%s\n\n" % (min_distance, "_"*140);		
					if(min_distance<maxdiff):
						models.append(mm);
						#closest = MarkovChain.merge(closest, mm);
					else:
						models.append(mm);
				else:
					models.append(mm)
					
				l = [];
			else:
				pass;
						
				
		pass
		
		
	def generate_string(self, length, chunk_size=None):
		'''Generates string on basis of the initiation and transition's probabilities in switching markov models'''
		seq = ''
		curlength = 0
		while(curlength<length):
			mm = weighted_choice_fast(self.models, self.switches)
			s = mm.generate_string()
			curlength += len(s)
			seq = "".join((seq, s));
			
			if(chunk_size and len(seq)>chunk_size):
				niter = len(seq)/chunk_size
				for i in range(niter):
					yield seq[chunk_size*i:chunk_size*(i+1)]
				seq = seq[chunk_size*niter:]
			
		terminal = len(seq) - (curlength-length)	
		yield seq[:terminal]
			
			
			
			
			

def markov_difference(model1, model2, pval=0.05):
	'''calculates the difference between two markov models'''
	emdiff = {};
	trdiff = {};

	for state in set(model1.emsupport.keys() + model2.emsupport.keys()):
		count1 = max(model1.emsupport.get(state, 1), 1)
		count2 = max(model2.emsupport.get(state, 1), 1)
		#print "%s\t%d\t%d\t%d\t%d" % (state, count1, count2, model1.support, model2.support)
		boundaries1 = wilson_score_interval_binomial(count1, model1.support, pval)
		boundaries2 = wilson_score_interval_binomial(count2, model2.support, pval)
		#print boundaries1, boundaries2
		s, e = numerictools.overlap(boundaries1, boundaries2)
		diff = max(0, s-e)
		emdiff[state] = diff**2;

	return sum(emdiff.values())/len(emdiff.values());			
			
			
			
if(__name__ == "__main__"):
	find_switch(1, 1, 1, 1, **conf)
	#seq = "ATATATTATATATTATATTATATATATATATATTTATTATTATFF"
	#order = 2;
	#mm = Markov_Model.from_sequence(seq, order)
	#print mm.generate_string(15)
	
	#for k, d in tr.iteritems():
		#print k, "\t", dict(d)
	#print
	#for k, (items, interval) in transitions.iteritems():
		#print k, "\t", items, "\t", interval 	