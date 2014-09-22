#! /usr/lib/python
'''Collections of filters for sam records'''
from sequencetools import entropy;


def chimera_left(ar, minright = 15):
	return (ar.rlen - ar.qend >= minright)
	
	
def repetitive(ar, min_entropy):
	pass
	# return if entropy(ar.seq[ar.qstart: ar.qend])