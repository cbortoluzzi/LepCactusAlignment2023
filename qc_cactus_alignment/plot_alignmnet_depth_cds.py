#!/usr/bin/env python


# Author : @cb46


import argparse
from ete3 import Tree
from pathlib import Path
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(description = 'Plot alignment depth of protein-coding sequence')
parser.add_argument('--d', help = 'Path to protein-coding alignment depth folders, one per species')
parser.add_argument('--tree', help = 'Phylogenetic tree')
parser.add_argument('--species_list', help = 'Tab delimited species list file')
parser.add_argument('--o', help = 'Output directory')



def order_species(tree):
	list_nodes = []
	t = Tree(tree)
	# Reroot the tree to the outgroup
	t.set_outgroup('Hydropsyche_tenuis')
	for node in t.traverse('postorder'):
		if node.is_leaf() and node.name != 'Hydropsyche_tenuis':
			list_nodes.append(node.name)
	return list_nodes


def get_species_group(species_list):
	superfamilies = {}
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, phylo_class, species_name, superfamily = line.strip().split()
			superfamilies[species_name] = superfamily
	return superfamilies



def alignment_depth_cds(list_tsv):
	mydict = {}
	for file in list_tsv:
		species_name = '_'.join(Path(file).stem.split('_')[0:2]).capitalize()
		genome = 0
		total_depth = 0
		with open(file) as f:
			for line in f:
				chromosome, start, end, length, depth = line.strip().split()
				length = int(length)
				depth = int(depth)
				genome += length
				total_depth += depth
		fraction_cds_aligned = round((total_depth / genome) * 100, 2)
		mydict[species_name] = fraction_cds_aligned
	return mydict


def plot_alignment_depth(list_nodes, superfamilies, mydict, output_directory):
	x, y = [], []
	for species in list_nodes:
		if species in mydict.keys():
			x.append(species)
			y.append(mydict[species])
		#else:
		#	x.append(species)
	# Manually create colors following Wright et al.
	list_colors = {'Noctuoidea': 'y', 'Bombycoidea': 'peru', 'Geometroidea': 'palevioletred', 'Drepanoidea': 'steelblue', 'Pyraloidea': 'gold', 'Papilionoidea': 'darkturquoise', 'Hesperioidea': 'darkgray',
	'Gelechioidea': 'coral', 'Zygaeinoidea': 'yellow', 'Cossoidea': 'slateblue', 'Torticoidea': 'yellowgreen', 'Tineoidea': 'cornflowerblue'}
	list_superfamily = [superfamilies[species] for species in x]
	colors = [list_colors[superfamily] for superfamily in list_superfamily]
	fig, ax = plt.subplots(figsize=(15, 8))
	plt.bar(x, y, color = colors, edgecolor = 'black')
	plt.xticks(rotation = 90, ha = 'right', fontsize = 8)
	plt.ylabel('Fraction of aligned coding sequence')
	plt.ylim(0, 100)
	figure = Path(output_directory, 'alignment_depth_CDS.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')




if __name__ == "__main__":
	args = parser.parse_args()
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	list_tsv = list(Path(args.d).rglob('*/*.tsv'))
	species_order = order_species(args.tree)
	species_superfamily = get_species_group(args.species_list)
	alignment_depth = alignment_depth_cds(list_tsv)
	plot_depth = plot_alignment_depth(species_order, species_superfamily, alignment_depth, args.o)

  
  
