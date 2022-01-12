#!/usr/bin/env python


# Author : @cb46


import argparse
from ete3 import Tree
from collections import defaultdict


parser = argparse.ArgumentParser(description = 'Prune phylogenetic tree')
parser.add_argument('-tree', help = 'Phylogenetic tree to prune')
parser.add_argument('-targetGenomes', help = 'List of species to be included in the pruned phylogenetic tree')
parser.add_argument('-pruned', help = 'Name of pruuned phylogenetic tree')



def prune_phylogenetic_tree(tree, target_genomes, pruned_tree):
	t = Tree(tree, format = 1)
	mydict = defaultdict(list)
	with open(target_genomes) as f:
		for line in f:
			species_name = line.strip().replace(' ', '_')
			for node in t.traverse('postorder'):
				if node.is_leaf():
					if node.name == species_name:
						mydict['prune'].append(node.name)

	for key in mydict:
		print (f"After pruning {len(mydict[key])} species remain")
		t.prune(mydict[key])
		t.write(format = 5, outfile = pruned_tree)



if __name__ == '__main__':
	args = parser.parse_args()
	pruned_phylo = prune_phylogenetic_tree(args.tree, args.targetGenomes, args.pruned)


	
