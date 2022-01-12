#!/usr/bin/env python


# Author : @cb46


import random
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description = 'Annotate phylogenetic tree')
parser.add_argument('--table', help = 'TSV file with information on class, clade, order, family, and group for each species, one per line')


def annotate_phylogenetic_tree(table):
	mydict = defaultdict(list)
	with open(table) as f, open('annotation_tree.txt', 'w') as out:
		next(f)
		for line in f:
			species, class_p, clade, order, family, group = line.strip().split('\t')
			species_name = species.replace(' ', '_')
			mydict[group].append(species_name)


		num_colors = len(mydict.keys())
		# Since in this case we only have two groups (i.e. butterfly and moth), we will select the color manually
		color = ['#ffffbf', '#998ec3']
		keys = list(zip(mydict.keys(), color))
		out.write('{}\n\n'.format('TREE_COLORS'))
		out.write('{}\n\n'.format('SEPARATOR SPACE'))
		out.write('{}\n\n'.format('DATA'))
		for i in keys:
			for j in mydict[i[0]]:
				key = ' '.join([j, 'range', i[1], i[0]])
				out.write('{}\n'.format(key))


if __name__ == "__main__":
	args = parser.parse_args()
	annotate_tree = annotate_phylogenetic_tree(args.table)
