'''this library is a collection of global constants(like path to the miRNA sequences, scripts and so on) used throughout the pipeline'''
import os;

root = "/data/BIO2/andrei/clash"
mroot = os.path.join(root, "mir")
sroot = os.path.join(root, "trunk")
proot = os.path.join(root, "projects")
rroot = os.path.join(root, "rscripts")
troot = os.path.join(root, "templates")

systems = ["hg19", "mm9", "ce6", 'v1']
#path to fasta files of miRNA
system2mir = {"hg19": os.path.join(mroot, "human.fa"), "mm9": os.path.join(mroot, "mouse.fa"), "ce6": os.path.join(mroot, "cel.fa"), "v1": os.path.join(mroot, "kshv_ebv_human.fa")}
#
experiment = {"hc": "HITS-CLIP", "pc": "PAR-CLIP"}

#bowtie options
bscore = 'C,30';
blength = 13;
bmistake = 0;

#conservation
conservation = {"hg19" : ["mm9", "hg19", "canFam2", "rn4"], 'ce6' : ['ce6', 'caeJap1', 'caePb2', 'caeRem3', 'cb3', 'priPac1']}

#html links
links = {"ce6": "<a href=\"http://141.80.186.52/cgi-bin/hgTracks?hgsid=7901&org=C.elegans&db=ce6&position=", "hg19": "<a href=\"http://141.80.186.52/cgi-bin/hgTracks?hgsid=7901&org=H.sapiens&db=hg19&position=", "mm9": "<a href=\"http://141.80.186.52/cgi-bin/hgTracks?hgsid=7901&org=M.musculus&db=mm9&position="}
