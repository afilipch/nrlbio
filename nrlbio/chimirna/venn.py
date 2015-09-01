from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles
import inspect 
import time
import argparse

parser = argparse.ArgumentParser(description='Script produce statistics and filtering of chiparts');
# input files
parser.add_argument('numbers', metavar = 'N', nargs = '+', type = str, help = "cardinalis of manifolds: only A, only B, only A&B, only C, only A&C, only B&C, A&C&B");
parser.add_argument('--labels', nargs = '?', type = str, help = "comma dived labels, should be enclosed by quotes");
args = parser.parse_args();

	
clash_numbers = (43, 119, 129, 197, 102, 219, 53)
clash_labels = ('ligation samples', 'control samples', 'normal samples')


plt.figure(figsize=(4,4), facecolor = "white")
if(args.labels):
	labels = args.labels.split(",");
	v = venn3(subsets=[int(x) for x in args.numbers], set_labels = labels, set_colors = ('0.40', "0.60", "0.80"))
else:
	v = venn3(subsets=[int(x) for x in args.numbers], set_colors = ('0.40', "0.60", "0.80"))
v.get_patch_by_id('111').set_color('0.05')
v.get_patch_by_id('011').set_color('0.45')
v.get_patch_by_id('101').set_color('0.25')
v.get_patch_by_id('110').set_color('0.35')
#v.get_patch_by_id('111').set_alpha(1.0)
#v.get_label_by_id('100').set_text('Unknown')
#v.get_label_by_id('A').set_text('Set "A"')
#c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
#c[0].set_lw(1.0)
#c[0].set_ls('dotted')
#plt.title("Sample Venn diagram")
#plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
            #ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.1, color = "white"))
plt.show(block = True)
#time.sleep(10) 
plt.close()