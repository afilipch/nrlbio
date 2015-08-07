#! /usr/lib/python
'''Creates makefile and directory structure for testing chiflex project'''
import argparse
import sys
import os

interaction_types = ['inter', 'intra', 'csj', 'lsj']

bowtie_settings = {'N': ('0','-'),
'L': ('16','-'),
'i': ('C,1', '-'),
'ignore-quals': ('True', '--'),
'norc': ('True', '--'),
'local': ('True', '--'), 
'mp': ('8,8', '--'),
'rfg': ('18,12', '--'),
'rdg': ('8,6', '--'),
'min-score': ('C,40', '--'), 
'k': ('4', '-'),
'D': ('30', '-'),
'R': ('2', '-'),
'no-unal': ('True', '--'), 
'f': ('False', '-'), 
'p': ('8', '-') }
bs_string = "\n".join(["\t%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in bowtie_settings.items()])

parser = argparse.ArgumentParser(description='Creates makefile and directory structure for testing chiflex projec');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the project folder. If folder does not exist, it will be created");
parser.add_argument('--reads', nargs = '?', type = str, required = True, help = "path to sequencing reads. fastq/fasta file");
parser.add_argument('--reference', nargs = '?', type = str, required = True, help = "path to the mapping reference");
parser.add_argument('--chiflex', nargs = '?', type = str, required = True, help = "path to the Chiflex folder")
parser.add_argument('--test', nargs = '?', type = str, required = True, help = "path to the test(folder with testing scrtipts for Chiflex) folder")
parser.add_argument('--name', nargs = '?', type = str, required = True, help = "name of the project, will be used as name for interactions ids")

parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "if set, a new makefile is created, but not folder structure");
parser.add_argument('--repetitive', nargs = '?', default = False, const=True, type = bool, help = "if set, repetitive mapped sequences are removed")
parser.add_argument('--nonunique', nargs = '?', default = False, const=True, type = bool, help = "if set, nonunique mappings with high alignment score are kept")

parser.add_argument('--bowtie', nargs = '+', default = [], type = str, help = "Bowtie settings. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones.\nDefault settings are:%s" % bs_string)
args = parser.parse_args();


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
	
bowtie_settings['x'] = os.path.abspath(os.path.abspath(args.reference)), '-'
bowtie_settings['U'] = os.path.abspath(os.path.abspath(args.reads)), '-'
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


def get_test_script(script, arguments={}, inp = '', out = None):
	'''Example: get_script(something.py, {'--a': ""}) will output: ('python [chiflex folder]/something.py'), '-a')
	'''
	l = [" ".join(('python', os.path.join(os.path.abspath(args.test), script) )), inp]
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
	input_files = os.path.abspath(args.reads)
	output_files = os.path.join('sam', 'mapped.sam')
	script = bs_list
	m.append(dependence(input_files, output_files, script))
	
	#Remove repetitive mappings if option 'repetitive' set
	if(args.repetitive):
		input_files = output_files
		output_files = os.path.join('sam', 'nonrep.bam')
		script = get_script('filter_sam.py', arguments={'--output': output_files, '--filters': 'repetitive', '--arguments': 'min_entropy=1.6'}, inp = input_files)
		m.append(dependence(input_files, output_files, script))
	
	#Collapse confident, but nonunique mappings into single sam entry, to protect them from further filtering out
	if(args.nonunique):
		input_files = output_files
		output_files = os.path.join('sam', 'collapsed.bam'), os.path.join('bed', 'collapsed.bed')
		script = get_script('collapse_nonunique_sam.py', arguments={'-s': output_files[0], '-b': output_files[1], '--minscore': 42}, inp = input_files)
		m.append(dependence(input_files, output_files, script))
		#following reassignment is done for furthe consistency
		output_files = output_files[0]
		
	#Demultiplex sam multiple hits into single hits, control single hits, chimeras and control chimeras
	input_files = output_files
	output_files = [os.path.join('sam', '%s.%s.bam' % (args.name, x)) for x in ['unique', 'unique_chimera', 'control_chimera', 'control']]
	script = get_script('demultiplex_chimera.py', arguments={'--output': 'sam', '--name': args.name, '--score': 'as_qstart', '--score_chimera': 'as_gap', '--maxgap': 8}, inp = input_files)
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
	script = get_script('filter_chimera.py', arguments={'-s': input_files[0], '-c' : input_files[1], '--features': 'AS1 AS2', '--fdr': 0.05}, out = output_files)
	m.append(dependence(input_files, output_files, script))
	
	#Filter single hits on basis of control single hits. LRG is applied for filtering
	input_files =  os.path.join('sam', '%s.unique.bam' % args.name) , os.path.join('sam', '%s.control.bam' % args.name) 
	output_files = os.path.join('sam', '%s.single_filtered.bam' % args.name) 
	script = get_script('filter_alignment.py', arguments={'-s': input_files[0], '-c' : input_files[1], '--features': 'AS qstart', '--fdr': 0.05, '--name': output_files})
	m.append(dependence(input_files, output_files, script))
	
	
	#Evalute mapping statistics
	input_files =  os.path.abspath(args.reads), os.path.join('sam', '%s.unique.bam' % args.name), os.path.join('sam', '%s.control.bam' % args.name), os.path.join('sam', '%s.unique_chimera.bam' % args.name) , os.path.join('sam', '%s.control_chimera.bam' % args.name), os.path.join('sam', '%s.nonunique.bam' % args.name), os.path.join('sam', '%s.nonunique_chimera.bam' % args.name), os.path.join('chimeras', 'filtered.bed'), os.path.join('sam', '%s.single_filtered.bam' % args.name)
	
	output_files = [os.path.join('evaluation_data', x)  for x in ['mapping_stat.yml', 'chimera_stat.yml', 'single_stat.yml', 'control_stat.yml', 'cf_stat.yml',  'sf_stat.yml']] 
	
	script = get_test_script('chiflex_evaluate.py', arguments={'-su': input_files[1], '-sc': input_files[2], '-cu': input_files[3], '-cc': input_files[4], '-sn': input_files[5], '-cn': input_files[6], '-cf' : input_files[7],'-sf': input_files[8]}, inp = input_files[0])
	m.append(dependence(input_files, output_files, script))
	
	
	#Visualize mapping statistics
	input_files = tuple(output_files)
	output_files = os.path.join('evaluation_plots', 'real2mapped', 'chimera.png')
	script = get_test_script('chiflex_visualize.py', arguments={'--mapping_stat': input_files[0], '--chimera_stat': input_files[1], '--single_stat': input_files[2], '--control_stat': input_files[3], '--fc_stat': input_files[4], '--fs_stat': input_files[5]})
	m.append(dependence(input_files, output_files, script))
	
	
			
	#makefile header
	if(isinstance(output_files, str)):
		m.insert(0, "SHELL=/bin/bash\n.DELETE_ON_ERROR:\n\nall: %s" % output_files );
	else:
		m.insert(0, "SHELL=/bin/bash\n.DELETE_ON_ERROR:\n\nall: %s" % " ".join(list(output_files)));
	# makefie cleaner
	m.append('clean:\n\techo "nothing to clean."\n');
	
	return "\n\n".join(m)
	
#create folder structure
project_path = os.path.abspath(args.path)
while (not args.only_makefile):
	try:
		os.makedirs(project_path);
		for folder in ('sam', 'chimeras', 'log'):
			os.mkdir(os.path.join(project_path, folder));
		break;	
	except:
		answer = raw_input("\nProject directory \'%s\' is currently exists, please type 'N' if you don't want to create a new project, 'MO' if you want to change/create only the Makefile, [another project name] if you want to create a new folder structure and makefile: " % project_path)
		if(answer=='N'):
			sys.exit('Project was not created')
		elif(answer=='MO'):
			sys.stderr.write('Makefile was changed/added to the existing project %s\n' % project_path)
			break
		else:
			project_path = os.path.abspath(answer)

	
with open(os.path.join(project_path, 'Makefile'), 'w') as mf:
	mf.write(makefile())
	
	
#report project call:
arguments_report = (
('name', ('Project name, assigned to the generated interactions', str)),
('path', ('Project folder', os.path.abspath)),
('reads', ('Sequencing reads', os.path.abspath)),
('chiflex', ('Chiflex module used in the project', os.path.abspath)), 
('reference', ('Reference(genome, transcriptome) used for mapping', os.path.abspath)), 
('repetitive', ('Repetitive mapped sequences are removed', str)),
('nonunique', ('Nonunique mappings with high alignment score are kept', str)), 
('only_makefile', ('New Makefile was generated', str)),  
)

with open(os.path.join(project_path, 'log/project.txt'), 'w') as rf:
	rf.write("Project call:\npython %s\n\n" % " ".join(sys.argv))
	
	for arg, (description, fun) in arguments_report: 
		av = getattr(args, arg);
		if(av):
			rf.write("%s:\t%s\n" % (description, fun(av)))
		else:
			rf.write("%s:\tnot set\n" % description)
			
	rf.write("bowtie2 settings:\n")
	for arg, (value, dashes) in bowtie_settings.items():
		rf.write("\t%s%s:\t%s\n" % (dashes, arg, value))
		
	