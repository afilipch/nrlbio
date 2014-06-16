# /usr/bin/python
'''generators for different kinds of input files. They may be used to allow multiprocessing and save memory footprint.
usually you should provide path to input file and size of elements(not lines) in buffer'''
import sys;
import itertools;
import sequencetools
from collections import *

fastq = namedtuple('fastq', 'id, seq, sign, qual')

def generator_fastq(path, take = ["id", "seq"], reverse = False, shuffle = False):
	'''yields fastq object(collections.namedtuple based) from fastq file provided
	
	path str: path to fastq file
	take list: list of fastq entry attributes to pass to yielded fastq object. Allows to save some memory. Entries have to be: "id", "seq", "sign", "qual"
	reverse bool: if True fastq sequence and quality string will be reversed
	shuffle bool: if True fastq sequence will be reversed
	
	Yields fastq: fastq object(collections.namedtuple based) 
	'''
	my_take = {"id": 0, "seq": 1, "sign": 2, "qual": 3}
	line_set = [];
	
	with open(path) as f:
		for line in f:  
			line_set.append(line.strip());
			if (len(line_set) == 4):
				appendix = [None]*4;
				for el in take:
					if(reverse and el in ["seq", "qual"]):
						appendix[my_take[el]] = line_set[my_take[el]][::-1];
					elif(shuffle and el == "seq"):
							appendix[my_take[el]] = sequencetools.shuffle_string(line_set[my_take[el]])
					else:	
						appendix[my_take[el]] = line_set[my_take[el]]
				yield fastq(*appendix)
				line_set = [];				

	
def generator_maf(path, aligned_species=None):
	from conservation import Maf;
	'''reads maf stitched file taking into account interactions id. Yields conservation.Maf object corresponding to each entry in maf file
	
	path string: path to maf stitched file
	aligned_species list: names of aligned species(genome suystems: ce6, mm9, etc.) to be considered in further analysis. If not given takes all species
	
	Yields conservation.Maf: object corresponding to each entry in maf file
	'''
	alignment = OrderedDict()
	with open(path) as f:
		for l in f:
			l = l.strip();  
			if(l):
				if(l[0] == ">"):
					if(len(l.split(":")) == 3):
						header = l
						specie = l.split(".")[0][1:];
					else:
						specie = l[1:];
				elif(aligned_species == None or specie in aligned_species):
					alignment[specie] = l.upper();
			else:
				yield Maf(header, alignment)
		try:
			yield Maf(header, alignment)
		except:
			pass

	
def grouper(iterable, n):
	it = iter(iterable);
	arr = [];
	while(it):
		for i in range(n):
			try:
				arr.append(next(it));
			except:
				yield arr;
				return
		yield arr;
		arr = []
				
#testing			
if(__name__ == "__main__"):
	#mygen = generator_fastq(sys.argv[1], take = ["id", "seq", "sign", "qual"], shuffle = True)
	#for fq in mygen:
		#print fq;
		#print;
	mygen = generator_maf(sys.argv[1])
	for fq in mygen:
		print fq;
		print;