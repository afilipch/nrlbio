# /usr/bin/python
'''collections of classes, wrappers and functions for html. Classes are designed to be called in applications to create objects, ready to pass in a Django template'''

class StatAttribute(object):
	'''class represents statistics regarding an attribute (age, length, score etc.) in a table like format
	Attributes:
		title str: name of attribue (age, length, score etc.)
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

def get_link(interval, system, internal = False):	
	if(internal):
		source = "http://141.80.186.52/"
	else:
		source = "http://genome.ucsc.edu/"	

	if(system == "ce6"):
		system_string = "cgi-bin/hgTracks?hgsid=7901&org=C.elegans&db=ce6&position=";
	elif(system == "hg19"):
		system_string = "cgi-bin/hgTracks?hgsid=7901&org=H.sapiens&db=hg19&position=";
	elif(system == "mm9"):
		system_string = "cgi-bin/hgTracks?hgsid=7901&org=M.musculus&db=mm9&position="; 
	else:
		pass;
		
	return "%s%s%s:%d..%d" % (source, system_string, interval.chrom , interval.start+1, interval.stop);
	
	