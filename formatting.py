# /usr/bin/python
'''collections of classes and functions to output python objects in a nice way'''

import collections;

def feature_dict_fraction(d, top = 0, key_names = None):
	'''ouputs formatted representation of feature_dict(any dictionary which keys are values of some feature)
	
	d dict: feature dictionary. Keys: values of certain feature, Values: support of corresponding key. For example {178: 25, 181: 13} may represent that case study where 25 people have height of 178cm and 13 of 181cm
	
	top int: how many entries of d to output. Entries are sorted according to the value of values. If top<1 outputs all entries.
	
	key_names iterable: in case of multifeatures (that is keys in d are iterables) helps to assign each name to each feature.
	
	Return: string. Representation of feature dictionary
	'''
	
	s = '';
	total = float(sum(d.values()))
	header = "fraction"	
	
	if(top < 1):
		top = len(d)
	else:
		pass;
			
	if(key_names):
		s += "\t".join([str(x) for x in key_names])
		s += "\t%s\n" % header;
	else:
		pass;
		
	for k, v in sorted(d.items(), key = lambda x: x[1], reverse = True)[:top]:
		vn = v/total;
		if (isinstance(k, collections.Iterable) and type(k) != str):
			s += "\t".join([str(x) for x in k])
		elif(k != ''):
			s += "%s\t" % str(k)
		else:
			s+= "%s\t" % str(None)
		s += "\t%1.5f\n" % vn
	return s

def feature_dict_total(d, top = 0, key_names = None):
	'''ouputs formatted representation of feature_dict(any dictionary which keys are values of some feature)
	
	d dict: feature dictionary. Keys: values of certain feature, Values: support of corresponding key. For example {178: 25, 181: 13} may represent that case study where 25 people have height of 178cm and 13 of 181cm
	
	top int: how many entries of d to output. Entries are sorted according to the value of values. If top<1 outputs all entries.
	
	key_names iterable: in case of multifeatures (that is keys in d are iterables) helps to assign each name to each feature
	
	Return: string. Representation of feature dictionary
	'''	
	s = '';
	header = "total"	
	
	if(top < 1):
		top = len(d)
	else:
		pass;
			
	if(key_names):
		s += "\t".join([str(x) for x in key_names])
		s += "\t%s\n" % header;
	else:
		pass;
		
	for k, v in sorted(d.items(), key = lambda x: x[1], reverse = True)[:top]:
		if (isinstance(k, collections.Iterable) and type(k) != str):
			s += "\t".join([str(x) for x in k])
		elif(k != ''):
			s += "%s\t" % str(k)
		else:
			s+= "%s\t" % str(None)
		s += "\t%d\n" % v
	return s	