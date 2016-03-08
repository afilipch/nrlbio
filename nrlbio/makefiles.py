'''Collection of functions supporting Makefile generation'''
import os
import sys

#chiflex_package = os.path.abspath(os.path.join(args.chiflex, 'chiflex'));
#splicing_package = os.path.abspath(os.path.join(args.chiflex, 'splicing'))
#interaction_package = os.path.abspath(os.path.join(args.chiflex, 'interaction'))
#mirna_package = os.path.abspath(os.path.join(args.chiflex, 'mirna'))

chiflex_package = r'~/nrlbio/chiflex'

def get_script(script, arguments={}, inp = '', out = None, package=chiflex_package):
	'''Example: get_script(something.py, {'--a': 7}, 'inp.txt', 'out.txt') will output: ('python [chiflex_package]/something.py'), 'inp.txt', '--a', '7', '>', 'out.txt')'''
	if(hasattr(inp, '__iter__')):
		input_files = " ".join(inp);
	else:
		input_files = inp
	
	l = [" ".join(('python', os.path.join(package, script) )), input_files]
	for k,v in arguments.items():
		if(v==False):
			pass;
		else:
			l.append(k)
			l.append(str(v))
	if(out):
		l.append(">")
		l.append(out)
	return l
		
	
def dependence(input_files, output_files, script):
	'''Creates Makefile dependence'''
	if(hasattr(input_files, '__iter__')):
		inp = " ".join(input_files);
	else:
		inp = input_files
	if(hasattr(output_files, '__iter__')):
		out = " ".join(output_files);
	else:
		out = output_files
	return "%s: %s\n\t%s" % (out,inp, " ".join([str(x) for x in script]))


def get_header(output_files, phony=False):
	if(isinstance(output_files, str)):
		ofs = output_files 
	else:
		ofs = " ".join(list(output_files))
	
	if(phony):
		return "SHELL=/bin/bash\n.DELETE_ON_ERROR:\n\nall: %s\n.PHONY: all %s" % (ofs, ofs)
	else:
		return "SHELL=/bin/bash\n.DELETE_ON_ERROR:\n\nall: %s" % ofs
	

def get_bowtie_call(settings, arguments, reference, reads, project_name):
	for bo in arguments:
		try:
			name, value = bo.split("=");
		except(ValueError):
			sys.exit("Bowtie options are provided in a malformatted way. Please see help for example\n") 
			
		if(name in settings):
			settings[name] = (value, settings[name][1])
		else:
			sys.stderr.write("provided option \'%s\' is currently not supported by Chiflex and will be ignored\n" % name) 
		
	settings['x'] = os.path.abspath(os.path.abspath(reference)), '-'
	settings['U'] = os.path.abspath(os.path.abspath(reads)), '-'
	settings['S'] = os.path.join('sam', '%s.mapped.sam' % project_name), '-'


	bs_list = ['bowtie2'];
	for k, v in settings.items():
		if(v[0]=='True'):
			bs_list.append(v[1] + k);
		elif(v[0]=='False'):
			pass;
		else:
			bs_list.append(v[1] + k);
			bs_list.append(v[0])
			
	return bs_list;


def get_bowtie_help(bowtie_configurations):
	bowtie_help_list = [];
	for mode, settings in bowtie_configurations.items():
		bowtie_help_list.append("\n%s:" % mode)
		bowtie_help_list.append("[%s]" % " ".join(["%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in settings.items()]))
	return " ".join(bowtie_help_list)	

