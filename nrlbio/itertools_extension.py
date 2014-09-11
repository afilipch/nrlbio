# /usr/bin/python
'''library supporting advanced handling, processing, conversion of iterables'''

import sys
import os
import itertools
from collections import defaultdict

class LengthException(Exception):
	pass;
	
	
class ArgumentException(Exception):
	pass;	

	
	
	
def pop_first(iterable):
	'''pops first element of any iterable, regardless does it has __getitem__ method or not'''
	if(hasattr(iterable, 'next')):	
		return iterable.next();
	else:
		return iterable.pop(0);

		
def iterable_of_objects_to_counts_dict(iterable, attributes=[]):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		attributes list: attributes used to pull objects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
	'''
	d = defaultdict(int);	
	
	if(hasattr(iterable, 'next')):	
		first = iterable.next();
		if(not attributes):
			attributes += vars(first).keys()			
		k = tuple([getattr(first, attr) for attr in attributes])
		d[k]+=1;	
	else:
		first = iterable[0];
		if(not attributes):
			attributes += vars(first).keys()	
	
	for el in iterable:
		k = tuple([getattr(el, attr) for attr in attributes])
		d[k]+=1;
	return d;
	
	
def iterable_of_lists_to_counts_dict(iterable, indices=[]):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes(corresponding to indices), Value: number of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		indices list: indices used to pull objects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes(corresponding to indices), Value: number of objects in iterable with values of attributes equal to the ones in Key.
	'''
	d = defaultdict(int);
	
	if(hasattr(iterable, 'next')):	
		first = iterable.next();
		length = len(first);		
		if(not indices):
			indices += range(length)			
		k = tuple([first[i] for i in indices])
		d[k]+=1;
	else:
		first = iterable[0];
		length = len(first);		
		if(not indices):
			indices += range(length)					
		
	for el in iterable:
		if(len(el)!=length):
			raise LengthException('elements in iterable are of different length\n')
		else:		
			k = tuple([el[i] for i in indices])
			d[k]+=1;
	return d;
	
	
def iterable_of_objects_to_list_dict(iterable, attributes=[]):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: list of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		attributes list: attributes used to pull objects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes, Value: list of objects in iterable with values of attributes equal to the ones in Key.
	'''
	d = defaultdict(list);

	if(hasattr(iterable, 'next')):	
		first = iterable.next();
		if(not attributes):
			attributes += vars(first).keys()			
		k = tuple([getattr(first, attr) for attr in attributes])
		d[k].append(first);	
	else:
		first = iterable[0];
		if(not attributes):
			attributes += vars(first).keys()	
			
	for el in iterable:
		k = tuple([getattr(el, attr) for attr in attributes])
		d[k].append(el);
	return d;
	
	
def iterable_of_lists_to_list_dict(iterable, indices=[]):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes(corresponding to indices), Value: list of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		indices list: indices used to pull objects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes(corresponding to indices), Value: list of objects in iterable with values of attributes equal to the ones in Key.
	'''	
	d = defaultdict(list);
	
	if(hasattr(iterable, 'next')):	
		first = iterable.next();
		length = len(first);		
		if(not indices):
			indices += range(length)			
		k = tuple([first[i] for i in indices])
		d[k].append(first);
	else:
		first = iterable[0];
		length = len(first);		
		if(not indices):
			indices += range(length)	
		
	for el in iterable:
		if(len(el)!=length):
			raise LengthException('elements in iterable are of different length\n')
		else:		
			k = tuple([el[i] for i in indices])
			d[k].append(el);
	return d;
	

def iterable_to_dict(iterable, entry, mode, attributes=[]):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects	
		
		entry list|object: type of iterable entry
			if list: entry will be treated as an iterable with __getitem__ method. Attributes has to be list of integers or None
			if object: entry will be treated as an object with __getattr__ method. Attributes has to be list of string or None
		
		mode list|count: 
			if list: entries with the same attributes will pulled in a list
			if count: number of entries with the same attributes will be stored
			
		attributes list: attributes used to pull objects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
		OR indices used to pull objects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes, Value: number or list of objects in iterable with values of attributes equal to the ones in Key.
	'''
	if(entry == 'object' and mode == 'count'):
		 return iterable_of_objects_to_counts_dict(iterable, attributes=attributes);
	elif(entry == 'object' and mode == 'list'):
		 return iterable_of_objects_to_list_dict(iterable, attributes=attributes);
	elif(entry == 'list' and mode == 'count'):
		 return iterable_of_lists_to_counts_dict(iterable, indices=attributes);
	elif(entry == 'list' and mode == 'list'):
		 return iterable_of_lists_to_list_dict(iterable, indices=attributes);
	else:
		raise ArgumentException('Wrong arguments were passed to \'iterable_to_dict\'\nentry has to set to \'list\' or \'object\'\nmode has to set to \'list\' or \'count\'\n')
		
	
	
#testing section
if(__name__ == "__main__"):
	from random import randint
	
	class A(object):
		def __init__(self):
			self.a = randint(1,10);
			self.b = randint(1, 5);
			self.c = bool(randint(0,1));
		
	it_objects = (A() for _ in range(1000))
	it_lists = [(randint(1,4), randint(1,4), randint(1,4)) for _ in range(1000)]
	
	for kv in iterable_of_objects_to_list_dict(it_objects, ['a', 'c']).iteritems():
		print "%s|\t\t%s\n" % (kv[0], "\t".join([str(vars(x)) for x in kv[1]]));
	
	for kv in iterable_of_lists_to_list_dict(it_lists, [0, 2]).iteritems():
		print "%s|\t\t%s\n" % (kv[0], "\t".join([str((x)) for x in kv[1]]));
		
	#it_lists.append((1,7));
	
	#for kv in iterable_of_lists_to_list_dict(it_lists, [0, 2]).iteritems():
		#print "%s|\t\t%s\n" % (kv[0], "\t".join([str((x)) for x in kv[1]]));	
	
	
	
	
	
	