# /usr/bin/python
'''collections of classes, wrappers and functions for html. Classes are designed to be called in applications to create objects, ready to pass in a Django template'''

import sys

import jinja2

'''set jinja2 enviroment'''
env = jinja2.Environment(loader=jinja2.PackageLoader('nrlbio', 'templates'))

class StatAttribute(object):
	'''class represents statistics regarding an attribute (age, length, score etc.) in a table like format
	Attributes:
		title str: name of attribute (age, length, score etc.)
		headers list: headers of collumns
		entries list: rows in the table. element is an entry of a cell
	'''	
	def __init__(self, title, headers, entries = []):
		self.title = title;
		self.headers = headers;
		self.entries = entries;
		
class Stat(object):
	'''class represents statistics regarding a set of attributes (age, length, score etc.) in a table like format
	Attributes:
		title str: name of statistics (e.g Statistics of citizens in Paris)
		attributes list: element represents statistics regarding an attribute (age, length, score etc.) in a table like format
	'''		
	def __init__(self, title, attributes = []):
		self.title = title;
		self.attributes = attributes
		
	@classmethod	
	def from_stat_object(cls, title, obj, attributes, headers=None, ordered=False, top_entries=20):
		'''Creates stat instance from attributes of given object
		
			title str: name of statistics (e.g Statistics of citizens in Paris)
			obj object: any kind of object with attributes representing some statistics
			attributes list of strings: each string is a name of an 'obj' attribute, which is itself dictionary representing some statistics
			headers list of strings: each string is a name used to be a header for statistics' table. The order and the length should correspond to the 'attributes' list. If None, 'attributes' will be used as 'headers'.
			ordered bool: if True, attributes are ordered according to number of instances(that is, if the most common age in a group of people, then it will be the first entry), if False, attributes are ordered according to the attribute value(that is, if the youngest in a group of people is 16, then it will be the first entry)
			top_entries int: number of top entries to store. For example, if there 10 persons of age 45; 8 of 43; 6 of 47 and top_entries=2. Then only the following info: 45->10, 43-> 8 will be stored, but not 47->6. This option is used to avoid long html tables. Is ignored when 'ordered'==False
		
		Creates stat instance from attributes of given objects
		'''
		if(not headers):
			headers = attributes;
		elif(len(headers) != len(attributes)):
			raise Exception('length of headers has to be equal to the length of attributes, or not provided at all\n')

		
		instance = cls(title)			
		for a, header in zip(attributes, headers):
			attribute = getattr(obj, a)
			if(attribute):
				html_attribute = StatAttribute(header, [a, "total number", "fraction"], entries = [])
				total = float(sum(attribute.values()))
				if(ordered):
					for k, v in sorted(attribute.items(), key = lambda x: x[1], reverse = True)[:top_entries]:
						f = v/total
						html_attribute.entries.append([k, "%d" % v, "%1.5f" % f]);
				else:
					for k, v in sorted(attribute.items(), key = lambda x: x[0], reverse = False):
						f = v/total
						html_attribute.entries.append([k, "%d" % v, "%1.5f" % f]);					
				instance.attributes.append(html_attribute);
			else:
				sys.stderr.write('attribute: %s does not exist or empty in the %s object\n' % (a, str(obj)))				
		return	instance;
		
			
	def report(self, output = None, template = "statistic_tables.html"):
		'''creates html representation of the Stat oblect
			
			output str: path to the ouptut fot html report. If not set, STDIN is used
			template str: name of the html template for statistics report
			
		creates html representation of the Stat oblect
		'''
		t = env.get_template(template);	
		if(output):
			with open(output, 'w') as f:
				f.write(t.render({"statistics": self}))
		else:		
			t.render({"statistics": self})
			
			

			
def get_link(interval, system, internal = False):
	if(internal):
		source = "http://genome.mdc-berlin.de/"
	else:
		source = "http://genome.ucsc.edu/"
		
	if(system in ("hg19", 'hg38')):
		system_string = "cgi-bin/hgTracks?hgsid=7901&org=H.sapiens&db=%s&position=" % system;
	elif(system in ("mm9", "mm10")):
		system_string = "cgi-bin/hgTracks?hgsid=7901&org=M.musculus&db=%s&position=" % system; 
	elif(system in ("ce6", "ce10", "ce11")):
		system_string = "cgi-bin/hgTracks?hgsid=7901&org=C.elegans&db=%s&position=" % system;
	elif(system == 'circ'):
		return 'http://circbase.org/cgi-bin/singlerecord.cgi?id=%s' % interval.chrom
		
	else:
		raise ValueError('Genome system \'%s\' is not supported\n' % system);
		
	return "%s%s%s:%d..%d" % (source, system_string, interval.chrom , interval.start+1, interval.stop);


def get_tarbase_link(mirid):
	return "http://mirtarbase.mbc.nctu.edu.tw/php/search.php?q=search_exact&searchword=%s" % mirid


	
	