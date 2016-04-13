#! /usr/bin/python
'''Downloads genomes listed in the provided html file''' 
import os
import sys
import re
import argparse
from subprocess import Popen,call


parser = argparse.ArgumentParser(description='Downloads genomes listed in the provided html file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the html file. File is an html table of genomes used in multi-alignment, can be copypasted from ucsc conseravtion track description");
parser.add_argument('--lookin', nargs = '+', default = [], type = str, help = "If set, script looks in the listed folders for genomes which are already listed and creates symlinks to them")
args = parser.parse_args();


def extract_html(fname):
	'''generates tuple: url link to the genome, name of the genome assembly'''
	with open(args.path) as fname:
		for line in fname:
			M = re.search(r'\<A HREF="(?P<url>\S+)"\s+TARGET\=_blank\>(?P<name>.*)\</A\>', line)
			if not M:
				print "No match",line
			else:
				#print M.group()
				yield M.group('name').split('/')[-1], M.group('url').replace('../','http://genome.ucsc.edu/')


def try_wget_queries(name,url):
	'''Tries different sources to download genome from'''
	yield 'wget -c http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.fa.gz' % (name,name)
	yield 'wget -c http://hgdownload-test.cse.ucsc.edu/goldenPath/%s/bigZips/%s.fa.gz' % (name,name)
	yield 'wget -c ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.2bit' % (name,name)
	yield 'wget -c ftp://hgdownload.cse.ucsc.edu/gbdb/%s/%s.2bit' % (name,name)
	
	
    
for name,url in extract_html('genomes.html'):

	gzipped = "%s.fa.gz" % name
	fasta = "%s.fa" % name
	twobit = "%s.2bit" % name

	if os.path.exists(gzipped):
		print ">>> skipping %s, already downloaded" % gzipped
		continue
	else:
		pass;
	
	if os.path.exists(fasta):
		print ">>> skipping %s, already downloaded" % fasta
		continue
	else:
		pass;	

	if os.path.exists(twobit):
		print ">>> skipping %s, already downloaded" % twobit
		continue
	else:
		pass;
	
	
	for base in args.lookin:
		gpath = os.path.join(base, gzipped)
		fpath = os.path.join(base, fasta)
		bpath = os.path.join(base, twobit)
		if(os.path.exists(gpath)):
			os.symlink(gpath, gzipped)
			print ">>> skipping %s, already downloaded to %s. Symlink is created" % (name, gpath)
			break;
		if(os.path.exists(fpath)):
			os.symlink(fpath, fasta)
			print ">>> skipping %s, already downloaded to %s. Symlink is created" % (name, fpath)
			break;		
		if(os.path.exists(bpath)):
			os.symlink(bpath, twobit)
			print ">>> skipping %s, already downloaded to %s. Symlink is created" % (name, bpath)
			break;
		
	else:
		print ">>>",name,url
		for cmd in try_wget_queries(name,url):
			print ">>> trying to retrieve %s, executing %s" % (name,cmd)
			res = call(cmd,shell=True)
			if res == 0:
				print ">>> SUCCESS!"
				break


    #if int(res) != 0:
        #print "error, trying alternative"
        #cmd = 'wget -c http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/chromFa.tar.gz -O - | tar -xzf - -O | pigz > %s' % (name,fasta)
        #print "retrieveing %s, executing %s" % (name,cmd)

        #res = call(cmd,shell=True)
        #if int(res) != 0:
            #os.remove(fasta)
            #print "error, still failed with",name
            

## use the HTML! ../ in URL equals genome.ucsc.edu