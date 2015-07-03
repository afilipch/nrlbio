# /usr/bin/python
'''Library of classes and functions related to Markov Models'''

import sys
import random
import itertools
from collections import defaultdict, deque
import copy

import numpy as np
import yaml

from nrlbio.random_extension import weighted_choice_fast, weighted2interval, list2interval
from nrlbio.statistics.functions import wilson_score_interval_binomial
from nrlbio import numerictools

CONF_MULTIMARKOV = {"step_size": 50000,
"slide_size": 1000,
'modeldiff': 0.012,
'maxlook': 100,
'jump': 300,
'jumpdiff': 0.01,
'jumplook': 400,
'verbose' : False
}

def find_switch(string, growing, sliding, curpos, slide_size=1200, modeldiff=0.01, maxlook=100, jump=300, jumpdiff=0.009, jumplook=400, verbose=False, **kwargs):
	'''Function tries to find an exact position of a swith between different markov models on a given string. For example, for "CCCCCCCCCCCCCCCCCCCCCCCC|TTTTTTTTTTTTTTTTTTTTTTTT" two different markov models should be generated with a switch exactly at "|". NOTE: Solution is merely suboptimal, since additional heuristics is used to speed up an algorithm. Output also depends on heuristic's parameters, so they shoud be either kept default or changed consciously.
	
		Input arguments:	
	
		string str: string to generate separate markov models on
		growing GrowingMarkovChain: markov model, which tries iteratively extend along the string. The switch between it and sliding model has to be found. Shoud be initialized prior call of this function, on a string left to the "string" provided as input.
		sliding SlidingMarkovChain: markov model, which slides "slide_size" chars ahead the growing one in order to calculate a difference between them stepwise. The switch between it and growing model has to be found. Shoud be initialized prior call of this function, on first "slide_size" chars of "string"
		curpos int: current position of a "string" on a global sequence. Global sequence iteratively switches between models, each switch corresponds to a single call of this function with a "string" being substring of a global sequence.
		
		Parameters:
		
		slide_size int: length of a sliding model. The bigger the length the faster the algorithm but less precise.
		modeldiff float: minimal difference between growing and sliding windows to stop search and output switch position and growing model. The bigger the modeldiff, the fewer models(but more different from each other) are generated.
		maxlook int: we don't wnat to find a switch at a first position with difference between models more than modeldiff. To prevent it we have to look forward to find the true global optimum(position with the largest difference). The bigger it is, the more precise the algorithm but also more slow.
		jump int: If the difference between growing and sliding model is reasonbly low, it makes sense not to iterate symbol by symbol, but jump "jump" symbols further. The bigger the jump, the faster the algorithm but less precise.
		jumpdiff float: If the difference between growing and sliding model is reasonbly low, it makes sense not to iterate symbol by symbol, but jump "jump" symbols further if difference is less than jumpdiff. The bigger the jumpdiff, the faster the algorithm but less precise.
		jumplook int: jump also can happend if there is no jump or model switch happens for [jumplook] number of iterations. The bigger the jumplook, the faster the algorithm but less precise
		verbose bool: if True, function outputs to STDERR information on workflow
		
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
		
		if(lookforward==maxlook and maxscore>modeldiff):
			growing.shrink(string[i-lookforward-growing.order*2:i])
			ngrowing = GrowingMarkovChain.from_string(string[i-lookforward: i+slide_size-lookforward], growing.order)
			nsliding = SlidingMarkovChain.from_string(string[i+slide_size-lookforward: i+slide_size*2-lookforward], growing.order)
			if(verbose):
				sys.stderr.write("switch to a new model happens at position (%d), cause the difference (%1.5f) is more than max difference %1.5f and lookforward is more or equal than %d\n" % 
			(i+curpos-lookforward, maxscore, modeldiff, lookforward))
			return growing, ngrowing, nsliding, i-lookforward+slide_size;
			
		elif(diff<jumpdiff or lookforward>jumplook):
			growing.grow_long(string[i: i+jump])
			sliding = SlidingMarkovChain.from_string(string[i+jump: i+slide_size+jump], sliding.order)
			if(verbose):
				if(diff<jumpdiff):
					sys.stderr.write("jump of length (%d) happens at position (%d), cause the difference (%1.5f) is less than jump difference %1.5f\n" % (jump, i+curpos, diff, jumpdiff))
				else:
					sys.stderr.write("jump of length (%d) happens at position (%d), cause the lookforward %d is more than jumplook %d\n" % (jump, i+curpos, lookforward, jumplook))
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
		
		
def substring2models(string, order, slide_size=1000, **kwargs):
	'''Iteratively generates Markov models on a given string
	
		string str: string to generate models on;
		slide_size int: length of a sliding model. The bigger the length the faster the algorithm but less precise. The bigger the jump, the faster the algorithm but less precise
		
	Returns list: list of local models generated on a given string
	'''
	models = []
	curpos=slide_size
	
	cgrowing = GrowingMarkovChain.from_string(string[:slide_size], order);
	csliding = SlidingMarkovChain.from_string(string[slide_size:slide_size*2], order);

	while(cgrowing):
		growing, cgrowing, csliding, pos = find_switch(string[curpos:], cgrowing, csliding, curpos, slide_size=slide_size, **kwargs)
		curpos += pos;
		if(growing):
			growing._to_serializable();
			models.append(growing);
	
	return models;
	
	
def seqrecord2models(seqrecord, order, step_size=50000, slide_size=1000, **kwargs):
	'''Iteratively generates Markov models on a given string
	
		seq_record Bio.SeqRecord: sequence record to generate models on;
		order int: defines the length of the state(for example markov chain with order of two will consider states of the sequence "AGAT" as "AG", "GA", "AT")
		step_size int: it is time efficient to split huge seq_records into pieces. 
		slide_size int: length of a sliding model. The bigger the length the faster the algorithm but less precise.
		
	Returns list: list of local models generated on a given seqrecord
	'''	
	
	models = []
	steps = len(seqrecord)/step_size-1;
	for step in range(steps):
		models += substring2models(str(seqrecord[step*step_size : (step+1)*step_size].seq.upper()), order, slide_size=slide_size, **kwargs);
	else:
		models += substring2models(str(seqrecord[steps*step_size:].seq.upper()), order, slide_size=slide_size, **kwargs);
		
	return models;
	




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
		else:
			self.length = None;
		
		self._fix_integrity()
		
		self.emsupport = defaultdict(int, dict([(x[0], sum(x[1].values())) for x in self.transitions.items()]))
		self.support = sum(self.emsupport.values())
		
		
	def _fix_integrity(self):
		'''All possible states should get transition and emission pseudocounts'''
		
		#define all possible states 
		self.states = set(self.initiation.keys() + self.transitions.keys());
		for d in self.transitions.values():
			self.states.update(d.keys());
		
		#check if we can transit from any state
		for state in self.states:
			if(state not in self.transitions.keys()):
				self.transitions[state] = defaultdict(int);
		
		#check if all possible transitions are allowed
		for fr, d in self.transitions.iteritems():
			for state in self.states:
				if(state not in d.keys()):
					self.transitions[fr][state]=1;
								
			for to, c in d.items():
				if(not c):
					self.transitions[fr][to]=1;
					
					
				
	def _to_serializable(self):
		tr = {};
		for fr, d in self.transitions.items():
			tr[fr] = dict(d);
			
		self.transitions = tr;
		#print tr.items()
		max_init = max(tr.items(), key=lambda x: sum(x[1].values()))[0]
		self.initiation = {max_init: 1};
		
		self.states = None;
		self.emsupport = None;
		self.support = None;
			
			
				
			
			
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
		
		
	def generate_string(self, chunk_size, length=None, leading_seq=''):
		'''Generates string on basis of the initiation and transition's probabilities'''
		if(length is None):
			if(self.length is not None):
				length = self.length
			else:
				raise ValueError("length attribute has to be set explicitly, cause models length is unknown\n\n")
			
		tr = {};	
		for k, d in self.transitions.iteritems():		
			tr[k] = weighted2interval(d.items());
		init = weighted2interval(self.initiation.items());
				
		last = weighted_choice_fast(*init);
		l = [leading_seq, last]
		curlength = self.order + len(leading_seq)
		genlength = -len(leading_seq);

		
		for _ in range(length/self.order):
			last = weighted_choice_fast(*tr[last]);
			l.append(last);
			curlength+=self.order
			if(curlength>=chunk_size and genlength+chunk_size<length):
				s = "".join(l);
				overhang = s[chunk_size:]
				l = [overhang]
				curlength = len(overhang);
				genlength += chunk_size;
				yield s[:chunk_size];
		else:
			s = "".join(l);
			nl = length - genlength;
			if(nl>0):
				yield s[:nl];
				
			#if(genlength +nl != self.length):
				#print 'bo'

				
				
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
		
	def _to_serializable(self):	
		super(GrowingMarkovChain, self)._to_serializable()
		self.last = None;
		
	def grow(self, rsymbol):
		#print self.transitions
		k = self.order;
		nlast = "".join((self.last[1:], rsymbol));
		
		self.emsupport[nlast[-k:]] += 1
		self.transitions[nlast[-k*2:-k]][nlast[-k:]] += 1
		self.support+=1
		self.length+=1
		
		self.last = nlast;
		#print self.transitions
		#print self.support
		
	def grow_long(self, seq):
		for rsymbol in seq:
			self.grow(rsymbol);
			
			
	def shrink(self, seq):
		#sys.exit();
		k = self.order;
		self.last = nlast = seq[:k*2];
		
		
		for s in seq[k*2:]:
			nlast = "".join((nlast[1:], s));
			self.emsupport[nlast[-k:]] -= 1
			self.transitions[nlast[-k*2:-k]][nlast[-k:]] -= 1
			self.support-=1
			self.length-=1
			
			
		
			
			
		
		
class MultiMarkov(object):
	'''Set of connected(through stochastic switches) Markov Models. May be useful for sequences with regions of vastly different transitions' probabilities.
	NOTE: all the models correspond to ONE sequence
	
	Attributes:
		models list of markov models: list of Markov Models, each representing specific region of the sequence. NOTE: all the models correspond to ONE sequence
		switches list of floats: list of probabilities(in interval format) for selection of the models
	'''
	
	
	def __init__(self, models, order, switches=None):
		self.models = models;
		self.order = order;
		if(switches):
			self.switches;
		else:
			self.switches = [x.length for x in models]

		
	@classmethod	
	def from_seqrecord(cls, seqrecord, order, **kwargs):
		'''creates multiple marcov models (MultiMarkov model) from given seq_record. Each model corresponds to specific local symbolic(nucleotide) composition on a given sequence.
		NOTE: The main aim of MultiMarkov modelling is to find different local specifity on a given sequence. For example, for "CCCCCCCCCCCCCCCCCCCCCCCC|TTTTTTTTTTTTTTTTTTTTTTTT" two different markov models should be generated with a switch exactly at "|". 
		NOTE: Solution is merely suboptimal, since additional heuristics is used to speed up an algorithm. Output also depends on heuristic's parameters, so they shoud be either kept default(conf dict in this library) or changed consciously.
		
			required arguments:
				seqrecord Bio.SeqRecord: sequence record to generate models on;
				order int: defines the length of the state(for example markov chain with order of two will consider states of the sequence "AGAT" as "AG", "GA", "AT")
				
			Optional arguments in **kwargs:
				step_size int: it is time efficient to split huge seq_records into pieces. default=50000
				slide_size int: length of a sliding model. The bigger the length the faster the algorithm but less precise. default=1000
				modeldiff float: minimal difference between growing and sliding windows to stop search and output switch position and growing model. The bigger the modeldiff, the fewer models(but more different from each other) are generated. default=0.01
				maxlook int: we don't wnat to find a switch at a first position with difference between models more than modeldiff. To prevent it we have to look forward to find the true global optimum(position with the largest difference). The bigger it is, the more precise the algorithm but also more slow. default=100
				jump int: If the difference between growing and sliding model is reasonbly low, it makes sense not to iterate symbol by symbol, but jump "jump" symbols further. The bigger the jump, the faster the algorithm but less precise. default=300
				jumpdiff float: If the difference between growing and sliding model is reasonbly low, it makes sense not to iterate symbol by symbol, but jump "jump" symbols further if difference is less than jumpdiff. The bigger the jumpdiff, the faster the algorithm but less precise. default=0.009
				jumplook int: jump also can happend if there is no jump or model switch happens for [jumplook] number of iterations. The bigger the jumplook, the faster the algorithm but less precise. default=400
				
		Returns: MultiMarkov model generated from a given seqrecord
		'''
 
		models = seqrecord2models(seqrecord, order, **kwargs)
		return cls(models, order)

		
	def serialize(self, output):
		obj = {'order': self.order, 'models': self.models, 'switches': self.switches}
		with open(output, 'w') as f:
			f.write(yaml.dump(obj, default_flow_style=False))
		return True;	
			
	@classmethod		
	def deserialize(cls, serialized):
		with open(serialized, 'r') as f:
			d = yaml.load(f)
			models, order = d['models'], d['order']
			for m in models:
				m._fix_integrity();
		return	cls(models, order)
	
		
		
	def generate_string(self, chunk_size, length=None):
		'''Generates string on basis of the initiation and transition's probabilities in switching markov models'''
		leading_seq=''
		if(length is None):
			for mc in random.sample(self.models, len(self.models)):
				change_leading=False
				for s in mc.generate_string(chunk_size, leading_seq=leading_seq):
					if(len(s)==chunk_size):
						yield s;
					else:
						leading_seq = s;
						change_leading = True;
				if(not change_leading):				
					leading_seq = ''		
			else:
				if(leading_seq):
					yield leading_seq;
					
					
		else:		
			interval = list2interval(self.switches)
			generated = 0;
			
			while(length>generated):
				mc = weighted_choice_fast(self.models, interval);
				change_leading=False
				for s in mc.generate_string(chunk_size, leading_seq=leading_seq):
					if(length<=generated+chunk_size):
						yield s[:length-generated]
						generated += chunk_size;
						break;
					elif(len(s)==chunk_size):
						yield s;
						generated += chunk_size;
					else:
						leading_seq = s;
						change_leading = True;
				if(not change_leading):	
					leading_seq = ''
			
			
			

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