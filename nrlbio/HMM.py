# /usr/bin/python
'''Library of classes and functions related to (Hidden) Markov Model'''

import sys
import random
from collections import defaultdict

from nrlbio.random_extension import weighted_choice_fast, weighted2interval

class MarkovException(Exception):
	pass;	


class Markov_Model(object):
	def __init__(self, transitions, emissions, order):
		self.transitions = transitions;
		self.emissions = emissions;
		self.order = order
		self.fix_integrity()
		
		
	@classmethod
	def from_sequence(cls, seq, order):
		if(len(seq)<=order*2):
			raise MarkovException("it is impossible to generate Markov Model of order %d from the sequence of length(%d) less then double order\n" % (order, len(seq)))
		else:
			tr = defaultdict(lambda: defaultdict(int));
			em = defaultdict(int)
			for i in range(len(seq)-2*order+1):
				tr[seq[i:i+order]][seq[i+order:i+2*order]]+=1
				em[seq[i:i+order]]+=1
				
			transitions = {};
			for k, d in tr.iteritems():
				transitions[k] = weighted2interval(d.items());
				
			emissions = weighted2interval(em.items())
				
			return cls(transitions, emissions, order);
			
			
	def fix_integrity(self):
		states = set()
		for k,(items, interval) in self.transitions.items():
			states.update(items);
		
		for k in states:
			if(k not in self.transitions.keys()):
				self.transitions[k] = self.emissions;
			
			
	def generate_string(self, length):
		seq = weighted_choice_fast(*self.emissions);
		last = seq;
		for _ in range(length/self.order):
			#print last
			#print self.transitions[last]
			last = weighted_choice_fast(*self.transitions[last])
			#print last
			#print 
			seq = "".join((seq, last));
		return seq[:length];
		
		
		
		
class Meta_Markov_Model(object):
	def __init__(self, models, order, window):
		self.models = models;
		self.order = order;
		self.window = window;
		
		
	@classmethod
	def from_sequence(cls, seq, order, window):
		models = []
		niter = max(len(seq)/window-1, 0);
		
		for i in range(niter):
			models.append(Markov_Model.from_sequence(seq[i*window:(i+1)*window], order));
		models.append(Markov_Model.from_sequence(seq[niter*window:], order))
		
		return cls(models, order, window);
		
		
	def generate_string(self, length):
		niter = max(length/self.window-1, 0);
		last_length = length - niter*self.window;
		for _ in range(niter):
			mm = random.choice(self.models);
			yield mm.generate_string(self.window);
		mm = random.choice(self.models);
		yield mm.generate_string(last_length);			
			
			
			
			
if(__name__ == "__main__"):
	seq = "ATATATTATATATTATATTATATATATATATATTTATTATTATFF"
	order = 2;
	mm = Markov_Model.from_sequence(seq, order)
	print mm.generate_string(15)
	
	#for k, d in tr.iteritems():
		#print k, "\t", dict(d)
	#print
	#for k, (items, interval) in transitions.iteritems():
		#print k, "\t", items, "\t", interval 	