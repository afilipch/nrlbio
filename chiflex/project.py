#! /usr/lib/python
'''Creates makefile and directory structure for chiflex project'''
import argparse
import sys
import os

bowtie_settings = {'N': ('0','-'),
'L': ('16','-'),
'i': ('C,2', '-'),
'ignore-quals': ('True', '--'),
'norc': ('True', '--'),
'local': ('True', '--'), 
'mp': ('6,6', '--'),
'rfg': ('12,12', '--'),
'rdg': ('6,4', '--'),
'min-score': ('C,40', '--'), 
'k': ('4', '-'),
'D': ('30', '-'),
'R': ('2', '-'),
'no-unal': ('True', '--'), 
'p': ('8', '-') }
bs_string = "\n".join(["\t%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in bowtie_settings.items()])

parser = argparse.ArgumentParser(description='assignes hybridization energy to interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the project folder");
parser.add_argument('--reads', nargs = '?', type = str, required = True, help = "path to collapsed reads fastq file");
parser.add_argument('--reference', nargs = '?', type = str, required = True, help = "path to the mapping reference");
parser.add_argument('--chiflex', nargs = '?', type = str, required = True, help = "path to the Chiflex folder")
parser.add_argument('--name', nargs = '?', type = str, required = True, help = "name of the project, will be used as name for interactions ids")

parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "if set, a new makefile is created, but not folder structure");
parser.add_argument('--reports', nargs = '?', default = False, const = True, type = bool, help = "if set, html reports will be produced");
parser.add_argument('--stratify', nargs = '?', default = False, const = True, type = bool, help = "if set, interactions will be stratified according to interaction type");

parser.add_argument('--annotation', nargs = '?', type = str, default = None, help = "path to an annotation file in gff format. If provided, interacting loci will be annotated");
parser.add_argument('--exons', nargs = '?', type = str, default = None, help = "path to a file of exonic regions in bed format. If provided, specific type will be assigned to each interaction");

parser.add_argument('--advanced', nargs = '+', choices = ['repetitive', 'nonunique'], default = [], type = str, help = "advanced processing steps to potentialy improve perfomance\n\trepetitive: removes low entropy(repetitive) sequences from analysis\n\tnonunique: keeps high quality, but nonunique mappings")
parser.add_argument('--bowtie', nargs = '+', default = [], type = str, help = "Bowtie settings. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones.\nDefault settings are:%s" % bs_string)
args = parser.parse_args();

#print sys.argv

#parse advanced options
nonunique = 'nonunique' in args.advanced
repetitive = 'repetitive' in args.advanced

#parse bowtie options
for bo in args.bowtie:
	try:
		name, value = bo.split("=");
	except(ValueError):
		sys.exit("Bowtie options are provided in a malformatted way. Please see help for example\n") 
		
	if(name in bowtie_settings):
		bowtie_settings[name] = (value, bowtie_settings[name][1])
	else:
		sys.stderr.write("provided option \'%s\' is currently not supported by Chiflex and will be ignored\n" % name) 
	
bowtie_settings['x'] = os.path.abspath(args.reference), '-'
bowtie_settings['U'] = os.path.abspath(args.reads), '-'
bowtie_settings['S'] = os.path.join('sam', 'mapped.sam'), '-'

bs_list = ['bowtie2'];
for k, v in bowtie_settings.items():
	if(v[0]=='True'):
		bs_list.append(v[1] + k);
	elif(v[0]=='False'):
		pass;
	else:
		bs_list.append(v[1] + k);
		bs_list.append(v[0])
#######################################################################################################################
		
def get_script(script, arguments={}, inp = '', out = None):
	'''Example: get_script(something.py, {'--a': ""}) will output: ('python [chiflex folder]/something.py'), '-a')
	'''
	l = [" ".join(('python', os.path.join(os.path.abspath(args.chiflex), script) )), inp]
	for k,v in arguments.items():
		l.append(k)
		l.append(str(v))
	if(out):
		l.append(">")
		l.append(out)
	return l;	
		
	
def dependence(input_files, output_files, script):
	'''creates Makefile dependence'''
	if(hasattr(input_files, '__iter__')):
		inp = " ".join(input_files);
	else:
		inp = input_files
	if(hasattr(output_files, '__iter__')):
		out = " ".join(output_files);
	else:
		out = output_files
	return "%s: %s\n\t%s" % (out,inp, " ".join([str(x) for x in script]))	
	
	
def makefile():
	m=[];
	
	#Map reads with bowtie2
	input_files = args.reads
	output_files = os.path.join('sam', 'mapped.sam')
	script = bs_list
	m.append(dependence(input_files, output_files, script))
	
	#Remove repetitive mappings if option 'repetitive' set
	if(repetitive):
		input_files = output_files
		output_files = os.path.join('sam', 'nonrep.bam')
		script = get_script('filter_sam.py', arguments={'--output': output_files, '--filters': 'repetitive', '--arguments': 'min_entropy=1.6'}, inp = input_files)
		m.append(dependence(input_files, output_files, script))
	
	#Collapse confident, but nonunique mappings into single sam entry, to protect them from further filtering out
	if(nonunique):
		input_files = output_files
		output_files = os.path.join('sam', 'collapsed.bam'), os.path.join('bed', 'collapsed.bed')
		script = get_script('collapse_nonunique_sam.py', arguments={'-s': output_files[0], '-b': output_files[1], '--minscore': 42}, inp = input_files)
		m.append(dependence(input_files, output_files, script))
		#following reassignment is done for furthe consistency
		output_files = output_files[0]
		
	#Demultiplex sam multiple hits into single hits, control single hits, chimeras and control chimeras
	input_files = output_files
	output_files = [os.path.join('sam', '%s.%s.bam' % (args.name, x)) for x in ['unique', 'unique_chimera', 'control_chimera', 'control']]
	script = get_script('demultiplex_chimera.py', arguments={'--output': 'sam', '--name': args.name, '--score': 'as_qstart', '--score_chimera': 'as_gap', '--maxgap': 6}, inp = input_files)
	m.append(dependence(input_files, output_files, script))

	#Merge sam hits into chimeras in doublebed format
	input_files = os.path.join('sam', '%s.unique_chimera.bam' % args.name) 
	output_files = os.path.join('chimeras', 'unique.bed') 
	script = get_script('merged2chimeras.py', arguments={}, inp = input_files, out = output_files)
	m.append(dependence(input_files, output_files, script))
	
	#Merge sam hits into chimeras in doublebed format for decoy(control) reference
	input_files = os.path.join('sam', '%s.control_chimera.bam' % args.name) 
	output_files = os.path.join('chimeras', 'control.bed') 
	script = get_script('merged2chimeras.py', arguments={}, inp = input_files, out = output_files)
	m.append(dependence(input_files, output_files, script))
	
	#Filter chimeras on basis of control chimeras. LRG is applied for filtering
	input_files =  os.path.join('chimeras', 'unique.bed') , os.path.join('chimeras', 'control.bed') 
	output_files = os.path.join('chimeras', 'filtered.bed') 
	script = get_script('filter_chimera.py', arguments={'-s': input_files[0], '-c' : input_files[1], '--features': 'AS1 AS2 gap', '--fdr': 0.05}, out = output_files)
	m.append(dependence(input_files, output_files, script))

	#Reassign chimeras coordinates from positions on genomic features(mapping reference) to genomic ones
	input_files = output_files
	output_files = os.path.join('chimeras', 'assigned.bed') 
	script = get_script('assign_coordinates.py', arguments={}, inp = input_files, out = output_files)
	m.append(dependence(input_files, output_files, script))

	#Merge chimeras into interactions.That is, sequenced chimeric reads which have two respectively overlaping parts are merged together
	input_files = output_files
	output_files = os.path.join('interactions', 'interactions.bed'), os.path.join('interactions', 'interactions2readid.bed')
	script = get_script('chimera2interaction.py', arguments={'-oi': output_files[0], '-od': output_files[1], '--name': args.name, '--distance': -16}, inp = input_files)
	m.append(dependence(input_files, output_files, script))
	output_files = output_files[0]
	
	#Annotate interacting loci
	if(args.annotation):
		input_files = output_files, os.path.abspath(args.annotation)
		output_files = os.path.join('interactions', 'annotated.gff')
		script = get_script('annotate_bed.py', arguments={'--annotation': input_files[1]}, inp = input_files[0], out=output_files)
		m.append(dependence(input_files, output_files, script))
		
		input_files = output_files
		output_files = os.path.join('interactions', 'interactions.ordered.bed')
		script = get_script('order.py', arguments={}, inp = input_files[0], out=output_files)
		m.append(dependence(input_files, output_files, script))


	#Annotate type of interaction
	if(args.exons):
		input_files = output_files, os.path.abspath(args.exons)
		output_files = os.path.join('interactions', 'interactions.annotated.gff')
		script = get_script('annotate_chimera.py', arguments={'--exons': input_files[1]}, inp = input_files[0], out=output_files)
		m.append(dependence(input_files, output_files, script));
	
	if(args.stratify):
		

	#makefile header
	m.insert(0, "SHELL=/bin/bash\n.DELETE_ON_ERROR:\n\nall: %s" % output_files);
	# makefie cleaner
	m.append('clean:\n\techo "nothing to clean."\n');
	
	return "\n\n".join(m)
	
#create folder structure
project_path = os.path.abspath(args.path)
while (not args.only_makefile):
	try:
		os.makedirs(project_path);
		for folder in ('mkdir', 'sam', 'chimeras', 'interactions', 'html', 'reports', 'bed', 'log'):
			os.mkdir(os.path.join(project_path, folder));
		break;	
	except:
		answer = raw_input("Project directory \'%s\' is currently exists, please type 'N' if you don't want to create a new project, 'MO' if you want to change/create only the Makefile, [another project name] if you want to create a new folder structure and makefile: " % project_path)
		if(answer=='N'):
			sys.exit('Project was not created')
		elif(answer=='MO'):
			sys.stderr.write('Makefile was changed/added to the existing project %s\n' % project_path)
			break
		else:
			project_path = os.path.abspath(answer)

	
with open(os.path.join(project_path, 'Makefile'), 'w') as mf:
	mf.write(makefile())
	
	