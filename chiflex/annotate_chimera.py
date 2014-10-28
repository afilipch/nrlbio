#! /usr/lib/python
'''assignes the type(linera/circular splice junctions, intra/inter-molecular interaction) to each chimera/interaction.''' 
import argparse
import sys;

import pysam;

from pybedtools import BedTool, Interval;



parser = argparse.ArgumentParser(description='assignes the type(linear/circular splice junctions, intra/inter-molecular interaction) to each chimera/interaction');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-e', '--exons', nargs = '?', default = "sam", type = str, help = "path to the bed file of transcripts\' exons");
args = parser.parse_args();

def list2interval(l):
	#sys.stderr.write("%s\n" % str(l))
	return Interval(l[0], int(l[1]), int(l[2]), name=l[3], score=l[4], strand=l[5])

def chimera2intervals(path, n):
	with open(path) as f:
		for l in f:
			a=l.strip().split("\t");
			yield list2interval(a[6*n:6*(n+1)])


class annotated_chimera(object):
	def __init__(self, left_interval, right_interval, left_intersection, right_intersection):
		self.left_interval = left_interval
		self.right_interval = right_interval
		self.left_intersection = left_intersection
		self.right_intersection = right_intersection
		
		
	def _linear_splice_junction(self, li, ri, distance):
		if(li.strand == '+'):
			if(int(ri.score) - int(li.score) == 1 and abs(self.left_interval.stop - li.stop)<=distance and abs(self.right_interval.start - ri.start)<=distance):
				return True;
			else:
				return False;
		else:
			if(int(li.score) - int(ri.score) == 1 and abs(self.left_interval.start - li.start)<=distance and abs(self.right_interval.stop - ri.stop)<=distance):
				return True;
			else:
				return False;
				
	def _circular_splice_junction(self, li, ri, distance):
		if(li.strand == '+'):
			if(self.left_interval.start>self.right_interval.stop and abs(self.left_interval.stop - li.stop)<=distance and abs(self.right_interval.start - ri.start)<=distance):
				return True;
			else:
				return False;
		else:
			if(self.right_interval.start>self.left_interval.stop and abs(self.left_interval.start - li.start)<=distance and abs(self.right_interval.stop - ri.stop)<=distance):
				return True;
			else:
				return False;		
				
	def annotate(self, linear_distance=4, circular_distance=4):
		rdict = dict([(x[3], x) for x in self.right_intersection])
		same_transcript = False;
		
		for li in self.left_intersection:
			ri = rdict.get(li[3], None);
			if(ri):
				same_transcript = True;
				if(self._linear_splice_junction(li, ri, linear_distance)):
					return 'lsj';
				if(self._circular_splice_junction(li, ri, circular_distance)):
					return 'csj';
			else:
				pass;
				
		if(same_transcript):
			return "intra"
		else:
			return "inter"
		
		


			

def generate_annotated_chimeras(intersection1, intersection2):
	iter_temp = iter(intersection2)
	
	curname = '';
	left_intersection = [];
	right_intersection = [];
	
	for i1 in intersection1:
		#print curname;
		if(i1.name == curname):
			left_intersection.append(list2interval(i1[6:12]));
			curleft = list2interval(i1[:6])
		elif(curname):
			i2 = iter_temp.next()
			while(i2.name == curname):
				right_intersection.append(list2interval(i2[6:12]))
				curright = list2interval(i2[:6])
				i2 = iter_temp.next()
			else:
				#sys.stderr.write("%s\n" % str(curleft))
				yield annotated_chimera(curleft, curright, left_intersection, right_intersection);
				curleft = list2interval(i1[:6])
				curright = list2interval(i2[:6])
				left_intersection = [list2interval(i1[6:12])];
				right_intersection = [list2interval(i2[6:12])];
				curname = curleft[3]
		else:
			curleft = list2interval(i1[:6])
			left_intersection = [list2interval(i1[6:12])];
			curname = curleft[3]
	else:
		yield annotated_chimera(curleft, curright, left_intersection, right_intersection);
			
first_intervals = BedTool(chimera2intervals(args.path, 0))
second_intervals = BedTool(chimera2intervals(args.path, 1))
exons = BedTool(args.exons)

intersection1 = first_intervals.intersect(exons, wao=True, s=True)
intersection2 = second_intervals.intersect(exons, wao=True, s=True)

#print intersection1
#print intersection2

chimeras = open(args.path)
for ac, ch in zip(generate_annotated_chimeras(intersection1, intersection2), chimeras):
	print "%s\t%s" % (ch.strip(), ac.annotate(linear_distance=10, circular_distance=0))
	
chimeras.close()	

























