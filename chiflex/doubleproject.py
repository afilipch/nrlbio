#! /usr/bin/python
'''Creates makefile and directory structure for chiflex double reference project'''
import argparse
import sys
import os

from nrlbio.config.config import load_config;
from nrlbio.makefiles import dependence, get_header, get_bowtie_call, get_script

conf = load_config('doublechiflex')




#Bowtie options preliminary parsing
bowtie_configurations_step1 = conf['bowtie']['first'];
bowtie_configurations_step2 = conf['bowtie']['secons'];

bowtie_help_list = [];
for mode, settings in bowtie_configurations.items():
	bowtie_help_list.append("\n%s:" % mode)
	bowtie_help_list.append("[%s]" % " ".join(["%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in settings.items()]))
bowtie_help_str = " ".join(bowtie_help_list)
	

parser = argparse.ArgumentParser(description='Creates makefile and directory structure for chiflex project')#, formatter_class = argparse.RawTextHelpFormatter);

#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the project folder. If folder does not exist, it will be created");
parser.add_argument('--reads', nargs = '?', type = str, required = True, help = "path to sequencing reads. fastq/fasta file");
parser.add_argument('--indices', nargs = 2, type = str, required = True, help = "path to the mapping reference bowtie2 indices");
parser.add_argument('--chiflex', nargs = '?', type = str, required = True, help = "path to the Chiflex folder")
parser.add_argument('--name', nargs = '?', type = str, required = True, help = "name of the project, will be used as name for interactions ids")


#Paths to the files for additional annotation;
parser.add_argument('--genome', nargs = '?', type = str, default = None, help = "path to the reference");
parser.add_argument('--exons', nargs = '?', type = str, default = None, help = "path to a file of exonic regions in bed format. If provided, specific type(circular or linear splice junction, intra- or intermolecular interaction) will be assigned to each interaction. That is, only known splice sites will be found. NOTE: if [--genome] and [--exons] are provided simultaneously, then splice sites will be de novo annotated based on [--genome]");
parser.add_argument('--annotation', nargs = '?', type = str, default = None, help = "path to an annotation file in gff format. If provided, found genomic loci will be annotated");
parser.add_argument('--genome_index', nargs = '?', type = str, default = None, help = "reference(genome) index file (.fai)");

#Options for the mapping result postprocessing
parser.add_argument('--repetitive', nargs = '?', default = False, const=True, type = bool, help = "if set, repetitive mapped sequences are removed")
parser.add_argument('--nonunique', nargs = '?', default = False, const=True, type = bool, help = "if set, nonunique mappings with high alignment score are kept")
parser.add_argument('--reassign', nargs = '?', default = False, const=True, type = bool, help = "if set, hits position on a reference will be reassigned to the genomic ones. Usefull in the case of nongenomic references(transcriptome, rRNAs, etc.). NOTE: reference headers has to be in [chrom]|[strand]|[start]|[stop] format")


#Options for the output control
#parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "if set, a new makefile is created, but not folder structure");
#parser.add_argument('--reports', nargs = '?', default = False, const = True, type = bool, help = "if set, html reports will be produced");

#bowtie2 options
parser.add_argument('--bowtie', nargs = '+', default = [], type = str, help = "Bowtie settings. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones. Default settings for the modes are:%s" % bowtie_help_str)
args = parser.parse_args();

#######################################################################################################################
#Set priotity for the modes
main_mode = min(args.modes, key= lambda x: modes_order[x]);


#######################################################################################################################
##Parse bowtie options
bowtie_settings = bowtie_configurations[main_mode];
bs_list = get_bowtie_call(bowtie_settings, args.bowtie, args.reference, args.reads, args.name)


#######################################################################################################################
#Configure input arguments
if('interaction' in args.modes or 'splicing' in args.modes):
	if(args.genome):
		annotate_with_genome = True;
	elif(args.exons):
		annotate_with_genome = False;
	else:
		raise argparse.ArgumentError(arg_modes, "Without setting [--genome] or [--exon] option interaction and splicing modes cannot be used\n")
		
if(('clustering' in args.modes) and (not args.genome_index)):
	raise argparse.ArgumentError(arg_modes, "Without setting [--genome_index] option clustering mode cannot be used\n")
	
	
chiflex_package = os.path.abspath(os.path.join(args.chiflex, 'chiflex'));
splicing_package = os.path.abspath(os.path.join(args.chiflex, 'splicing'))
interaction_package = os.path.abspath(os.path.join(args.chiflex, 'interaction'))