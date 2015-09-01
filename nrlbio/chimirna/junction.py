import parsing;
import sys;
import os;
import argparse;
import chipart_lib;
import gconst;
from collections import *;

parser = argparse.ArgumentParser(description='reanalize partition of the reads into miRNA and target parts');
#parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to left chipart file (adjusted_chiparts.tsv)");
parser.add_argument('-r', '--right', nargs = '+', type = str, required = True, help = "path file of filtered right chiparts(sam/right.tsv)");
parser.add_argument('-l', '--left', nargs = '+', type = str, required = True, help = "path to left chipart file (adjusted_chiparts.tsv)");
parser.add_argument('-f', '--fastq', nargs = '+', type = str, help = "path to potential chimeras in fastq format(left/filtered_candidates.fastq)");
parser.add_argument('-g', '--gmode', nargs = '?', default = False, const = True, help = "account fot prefernce of T1 RNAse to cut after G in PAR-CLIP protocol");
parser.add_argument('-c', '--tc', nargs = '?', default = False, const = True, help = "allows t->c conversions in target extension");
parser.add_argument('-d', '--depth', nargs = '?', default = 10, type = int, help = "depth of search");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
args = parser.parse_args();

system = args.system
exec("from sequence_data.systems import %s as gsys" % system);

Right = namedtuple("Right", 'chromosome, genome_start, genome_end, read_id, score, strand, read_start, read_end, collapse, conversion, aligned')

mirdict = parsing.fasta2dict(gconst.system2mir[system]);

def adjust_tc(rp):
	conv1 = rp.conversion - rp.read_start;
	if(conv1 < 0 or conv1 >= len(rp.aligned) - rp.aligned.count("-")):
		return rp.aligned, True;
	stop = 0;
	conv2 = 0;
	for el in rp.aligned:
		if(stop == conv1):
			break;
		else:	
			conv2 += 1;
			if(el != "-"):
				stop += 1
			else:
				pass;
	return rp.aligned[:conv2] + "C" + rp.aligned[conv2+1:], False;
			
	

def reanalize(left, right, seq, gmode, tc, depth):
	
	
	
	if(right.strand == "+"):
		refseq = gsys.genome.get_oriented(right.chromosome, right.genome_start-depth, right.genome_end, right.strand).upper()
	else:
		refseq = gsys.genome.get_oriented(right.chromosome, right.genome_start, right.genome_end+depth, right.strand).upper()	
		
	if(tc):
		aligned, notc = adjust_tc(right)
	else:
		aligned = right.aligned	
	
	#if(left.read_id == "@ago2_kishore-pc_9202443_x1"):
		#lseq = seq[left.read_start: left.read_end + right.read_start]
		#print left;
		#print refseq		
		#print lseq;
		#extension = 0;
		#gextension = -5;
		#for i in range(depth):
			#n1 = lseq[-1-i];
			#n2 = refseq[depth - i - 1];
			#print n1, n2
			#if(n1 == n2):
				#if(n1 == "G"):
					#gextension = i;
			#elif(tc  and notc and n2 == "T" and n1 == "C"):
				#notc = False;
				#print "b"
			#else:
				#extension = i;
				#break;	
		#else:
			#extension = i;	
		#print extension, gextension	
			
	
	
	lseq = seq[left.read_start: left.read_end]
	extension = 0;
	gextension = -5;
	#print lseq, left.read_id
	for i in range(depth):
		n1 = lseq[-1-i];
		n2 = refseq[depth - i - 1];
		if(n1 == n2):
			if(n1 == "G"):
				gextension = i;
		elif(tc  and notc and n2 == "T" and n1 == "C"):
			notc = False;
		else:
			extension = i;
			break;	
	else:
		extension = i;	
		
	if(gmode and gextension + 4 > extension):
		extension = gextension;
		
	#if(right.read_start != 0):
		#extension = 0;
		
			
	if(right.strand == "+"):
		nr = Right(right.chromosome, right.genome_start - extension, right.genome_end, right.read_id, right.score, right.strand, right.read_start, right.read_end + extension, right.collapse, right.conversion, lseq[len(lseq)-extension:] + aligned)
	else:
		nr = Right(right.chromosome, right.genome_start, right.genome_end + extension, right.read_id, right.score, right.strand, right.read_start, right.read_end + extension, right.collapse, right.conversion, lseq[len(lseq)-extension:] + aligned)	
		
	nl = chipart_lib.Chipart(left.read_id, left.ref_id, left.read_start, left.read_end-extension, left.read_length, left.ref_start, left.ref_end - extension, left.ref_length, left.cut_left, lseq[-1-extension], left.conversions, left.weight, left.score)	

	return nl, nr, extension;
		
rid2seq = {};
for p in args.fastq:
	f = open(p);
	lset = [];
	for l in f:
		lset.append(l.strip());
		if(len(lset) == 4):
			rid2seq[lset[0]] = lset[1];
			lset = [];
	f.close();

rid2rparts = {};
for p in args.right:
	f = open(p);
	for l in f:
		a = l.strip().split("\t");
		rid2rparts[a[3]] = Right(a[0],int(a[1]), int(a[2]), a[3], int(a[4]), a[5], int(a[6]), int(a[7]), int(a[8]), int(a[9]), a[10]);
	f.close();

chiparts = chipart_lib.read(args.left); 
rid2lpart = {};
for ch in chiparts:
	rid2lpart[ch.read_id] = ch;
	
fh_nl = open("left.tsv", 'w')	
fh_nr = open("right.tsv", 'w')	

cc = defaultdict(int);
for rid, right in rid2rparts.iteritems():
	left = rid2lpart[rid]
	seq = rid2seq[rid]
	if(left.read_end-left.read_start< args.depth):
		continue;
	if(right.read_start == 0):
		nl, nr, e = reanalize(left, right, seq, args.gmode, args.tc, args.depth);
		cc[e] += 1;	
	else:
		nl = left
		if(args.tc):
			aligned, notc = adjust_tc(right)
		else:
			aligned = right.aligned	
		nr = Right(right.chromosome, right.genome_start, right.genome_end, right.read_id, right.score, right.strand, right.read_start, right.read_end, right.collapse, right.conversion, aligned)	
		
	rarr = nr.chromosome, nr.genome_start, nr.genome_end, nr.read_id, nr.score, nr.strand, nr.read_start, nr.read_end, nr.collapse, nr.conversion, nr.aligned
	fh_nr.write("\t".join([str(x) for x in rarr]) + "\n");	
	
	larr = nl.read_id, nl.ref_id, nl.read_start, nl.read_end, nl.read_length, nl.ref_start, nl.ref_end, nl.ref_length, nl.cut_left, nl.cut_right, "|".join(nl.conversions), nl.weight, nl.score
	fh_nl.write("\t".join([str(x) for x in larr]) + "\n");		
	
fh_nl.close()
fh_nr.close()


#>>>>>>>>>>>>>>>>>>>>>>>> statistics output

ml = open(os.path.join("final", "mir_length.tsv"), 'w')
mc = open(os.path.join("final", "cut_length.tsv"), 'w')

r = Counter([x[1].ref_end - x[1].ref_start for x in rid2lpart.items()])
for i in range(min(r.keys())+1, max(r.keys())+1):
	ml.write("%d\t%1.4f\n" % (i, float(r[i])/sum(r.values())))

r = Counter([x[1].ref_length - x[1].ref_end for x in rid2lpart.items()])
for i in range(min(r.keys()), max(r.keys())+1):
	mc.write("%d\t%1.4f\n" % (i, float(r[i])/sum(r.values())))	

ml.close()
mc.close()





	