import os;
from gconst import *;

join = os.path.join;

def dependence(inp, out, arr):
	if(isinstance(out, list)):
		out = " ".join(out)
	if(isinstance(inp, list)):
		inp = " ".join(inp)		
	return "\n%s: %s\n\t" % (out,inp) + " ".join([str(x) for x in arr])

class Project(object):
	def __init__(self, name, reads, ref, index, system, mode, gmode, explore):
		self.name, self.reads, self.ref, self.index, self.system, self.mode, self.gmode, self.explore = name, os.path.abspath(reads), os.path.abspath(ref), os.path.abspath(index), system, mode, gmode, explore;
		self.mir = join(mroot, system2mir[system])
		if(explore):
			self.amode = "exp";
		else:
			self.amode = "hc";
		
		if(gmode):
			self.ginv = "--gmode"
		else:
			self.ginv = ""
			
		if(mode == "hc"):
			self.ft = "!!";
			self.crosslink = ""
		elif(mode == "pc"):
			self.ft = "CT"
			self.crosslink = 'plots/crosslink.pdf'
	def info(self):
		return "name of project: %s\npath to initial reads: %s\npath to mapping reference bed file: %s\npath to mapping reference index: %s\nmodel system: %s\ntype of experiment: %s\ngmode enabled: %r\nexplore mode enabled: %r\n" % (self.name, self.reads, self.ref, self.index, self.system, experiment[self.mode], self.gmode, self.explore)
	def makefile(self):
		r = 'SHELL=/bin/bash\n.DELETE_ON_ERROR:\n\nall: output/interactions.bed rstatistics/modes.tsv rstatistics/energy.tsv rstatistics/family.tsv plots/modes.pdf plots/energy.pdf plots/family.pdf plots/pattern.pdf plots/conversion.pdf %s'		% self.crosslink

		i = self.reads;
		o = 'left/anchors.tsv'
		l = ['python', join(sroot,'anchor.py'), i, '--ref', self.mir, '--mode', self.amode]
		r += dependence(i, o, l);
		
		i = self.reads;
		o = 'left/anchors_control.tsv';
		l = ['python', join(sroot,'anchor.py'), i, '--ref', self.mir, '--mode', self.amode, '--control']
		r += dependence(i, o, l);
		
		i = 'left/anchors.tsv'
		o = 'left/long.tsv'
		l = ['python', join(sroot,'length_statification.py'), i]
		r += dependence(i, o, l);
		
		i = 'left/anchors_control.tsv'
		o = 'left/long_control.tsv'
		l = ['python', join(sroot,'length_statification.py'), i, '-p']
		r += dependence(i, o, l);		
		
		i = ['left/long.tsv', 'left/long_control.tsv']
		o = 'right/right.fastq';
		l = ['python', join(sroot,'chipart.py'),  i[0],  '--control', i[1], '--fastq',  'left/candidates.fastq', self.ginv]
		r += dependence(i, o, l);
		
		i = 'right/right.fastq'
		o = 'right/control.fastq';
		l = ['python', join(sroot,'permute.py'),  i,  '-c', 8, '-d',  '>', o]
		r += dependence(i, o, l);
		
		i = 'right/right.fastq'
		o = 'right/extended.fastq';
		l = ['python', join(sroot,'introduce_conversion.py'),  i,  '-f', 'qfa', '--From', self.ft[0], '--to', self.ft[1],  '>', o]
		r += dependence(i, o, l)	
		
		i = 'right/control.fastq'
		o = 'right/extended_control.fastq';
		l = ['python', join(sroot,'introduce_conversion.py'),  i,  '-f', 'qfa', '--From', self.ft[0], '--to', self.ft[1],  '>', o]
		r += dependence(i, o, l)	
		
		i = 'right/extended.fastq'
		o = 'sam/true.sam';
		l = ['bowtie2', '-x', self.index, '-U', i, '--min-score', bscore, '-S', o, '-N', bmistake, '-L', blength, '-p', 5, '-i', 'C,1', '-D', 60, '-R', 4, '--local', '--reorder', '--no-unal', '--no-hd', '--no-sq', '--norc', '-k', 3]  
		r += dependence(i, o, l)	
		
		i = 'right/extended_control.fastq'
		o = 'sam/control.sam';
		l = ['bowtie2', '-x', self.index, '-U', i, '--min-score', bscore, '-S', o, '-N', bmistake, '-L', blength, '-p', 5, '-i', 'C,1', '-D', 60, '-R', 4, '--local', '--reorder', '--no-unal', '--no-hd', '--no-sq', '--norc', '-k', 3]  
		r += dependence(i, o, l)

		i = ['sam/true.sam', 'sam/control.sam']
		o = 'sam/mapreads.tsv';
		l = ['python', join(sroot,'mapreads.py'),  i[0],  '-c', i[1], '--ref', self.ref, '-s', self.system]
		r += dependence(i, o, l)	
		
		i = 'sam/mapreads.tsv';
		o = 'output/interactions.bed'
		l = ['python', join(sroot,'clustering.py'),  i,  '--chiparts',  'left/filtered_chiparts.tsv', '--mir', self.mir, '--fam', 'left/fam_ids.tsv', '-s', self.system]
		r += dependence(i, o, l)
		
		i = 'output/interactions.bed'
		o = 'rstatistics/modes.tsv'
		l = ['python', join(sroot,'modes.py'),  i]
		r += dependence(i, o, l)
		
		i = 'output/interactions.bed';
		o = 'rstatistics/energy.tsv';
		l = ['python', join(sroot,'hybrid.py'),  i];
		r += dependence(i, o, l);

		i = 'output/interactions.bed';
		o = 'rstatistics/family.tsv';
		l = ['python', join(sroot,'family.py'),  i, '-s', self.system];
		r += dependence(i, o, l);
		
		i = 'rstatistics/family.tsv';
		o = 'plots/family.pdf';
		l = ['Rscript', join(rroot,'family.r'),  i, o];
		r += dependence(i, o, l);
		
		i = 'rstatistics/energy.tsv';
		o = 'plots/energy.pdf';
		l = ['Rscript', join(rroot,'energy.r'),  i, o];
		r += dependence(i, o, l);	
		
		i = 'rstatistics/pattern.tsv';
		o = 'plots/pattern.pdf';
		l = ['Rscript', join(rroot,'pattern.r'),  i, o];
		r += dependence(i, o, l);
		
		i = 'rstatistics/modes.tsv';
		o = 'plots/modes.pdf';
		l = ['Rscript', join(rroot,'modes.r'),  i, o];
		r += dependence(i, o, l);
		
		i = 'rstatistics/conv_types.tsv';
		o = 'plots/conversion.pdf';
		l = ['Rscript', join(rroot,'conversion.r'),  i, o];
		r += dependence(i, o, l);
		
		if(self.mode == "pc"):
			i = 'rstatistics/crosslink.tsv';
			o = 'plots/crosslink.pdf';
			l = ['Rscript', join(rroot,'crosslink.r'),  i, o];
			r += dependence(i, o, l);		
		
		r += '\n\nclean:\n\techo "nothing to clean."\n'
		
		return r;