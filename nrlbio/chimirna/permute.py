#! /usr/bin/python
import random;
import argparse;
import sys;

parser = argparse.ArgumentParser(description='permute sequences in fa or qfa');
parser.add_argument('path', metavar = 'N', nargs = 1, type = str, help = "Path to reads");
parser.add_argument('-d', '--dinucleotide', nargs = '?', type = bool, default = False, const = True, help = "using this option preserve dinucleotide composition")
parser.add_argument('-c', '--coincide', nargs = '?', type = int, default = 0, const = 10, help = "using this option we disallow  x mer concidance between initial and permuted read")
args = parser.parse_args();

def get_seq_qual(init_seq, init_qual, dinucleotide):
	seq = ""
	qual = ""
	if(dinucleotide):
		a = range(0, (len(init_seq)+1)/2);
		random.shuffle(a);
		for i in a:
			seq += init_seq[i*2: i*2+2];
			qual += init_qual[i*2: i*2+2];  
		return(seq, qual);
	else:
		a = range(0, len(init_seq));
		random.shuffle(a);
		for i in a:
			seq += init_seq[i];
			qual += init_qual[i];  
		return(seq, qual);

def wrap_get_seq_qual(init_seq, init_qual, dinucleotide, coincide):
	count = 0;
	if(not coincide):
		return get_seq_qual(init_seq, init_qual, dinucleotide);
	else:
		my_coincide = True;
		trials = 0;
		while(my_coincide and trials < 4):
			ans = get_seq_qual(init_seq, init_qual, dinucleotide);
			my_coincide =  False;
			for i in range(0, len(init_seq) - coincide + 1):
				piece = init_seq[i: coincide + i]
				if(ans[0].find(piece) > -1):
					my_coincide = True;
					trials += 1;
				if (trials == 3):
					count = 1;
					break;
	return (ans, count);



count1 = 0;
count2 = 0;
fastq = file(args.path[0], 'r');  
line_set = [];
for line in fastq:  
	line = line.strip();
	line_set.append(line);
	if (len(line_set) == 4):
		ans = wrap_get_seq_qual(line_set[1], line_set[3], args.dinucleotide, args.coincide);
		seq, qual = ans[0]
		count1 += ans[1];
		count2 += 1;
		print line_set[0];
		print seq;
		print line_set[2]
		print qual;	
		line_set = [];          
fastq.close();
sys.stderr.write("total unpermuted:\t" + str(count1) + "\n")
sys.stderr.write("total reads:\t" + str(count2) + "\n")
  



