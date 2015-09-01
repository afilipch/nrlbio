import sys;
import itertools;
import miscellaneous;
#generators for different kinds of input files. They may be used to allow multiprocessing and save memory footprint.
#usually you should provide path to input file and size of elements(not lines) in buffer;


#>>>>>>>>>>>>> 
def generator_fastq(paths, take = ["id", "seq"], reverse = False, shuffle = False):
	c = 0;
	my_take = {"id": 0, "seq": 1, "sign": 2, "qual": 3}
	for path in paths:
		handler = open(path, 'r');  
		line_set = [];
		for line in handler:  
			line = line.strip();
			line_set.append(line);
			if (len(line_set) == 4):
				c += 1;
				appendix = [];
				for el in take:
					if(reverse and el in ["seq", "qual"]):
						appendix.append(line_set[my_take[el]][::-1])
					elif(shuffle and el == "seq"):
						appendix.append(miscellaneous.shuffle_str(line_set[my_take[el]]))						
					else:	
						appendix.append(line_set[my_take[el]])	
				yield tuple(appendix)
				line_set = [];
		handler.close(); 
	sys.stderr.write("total reads read %d\n" % c)

	
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
				
			

