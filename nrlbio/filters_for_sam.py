#! /usr/lib/python
'''Collections of filters for sam records'''
from sequencetools import entropy, chunk_entropy;


def chimera_left(ar, minright = 15):
	return (ar.rlen - ar.qend >= minright)
	
	
def repetitive(ar, min_entropy=1.4):
	return entropy(ar.query) > min_entropy
	
def contain_repetitive(ar, min_entropy=1.4):
	return chunk_entropy(ar.query, length=14, step = 2, order = 1) > min_entropy	
	