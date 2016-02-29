#! /usr/lib/python
'''Collections of filters for bed/gff/doublebed records'''
from nrlbio.numerictools import overlap


def distance(doublebed, maxoverlap = 1):
	mo = int(maxoverlap)
	start, end = overlap((doublebed[0].start, doublebed[0].stop), (doublebed[1].start, doublebed[1].stop))
	return end-start<=mo
	
	
def itype(doublebed, type_ = None):
	return (not type_) or (doublebed[0].attrs['interaction'] == type_)
	