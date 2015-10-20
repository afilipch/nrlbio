#! /usr/lib/python
'''Collections of filters for sam records'''
from nrlbio.sequencetools import entropy, chunk_entropy;


def chimera_left(ar, minright = 15):
	mr = float(minright)
	return (ar.rlen - ar.qend >= mr)
	
	
def repetitive(ar, min_entropy=1.4):
	me = float(min_entropy)
	return entropy(ar.query) > me
	
def contain_repetitive(ar, min_entropy=1.5, length=18):
	me = float(min_entropy)
	return chunk_entropy(ar.query, length, step = 2, order = 1) > me	
	