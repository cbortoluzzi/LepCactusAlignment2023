#!/usr/bin/env python


# Author : @cb46


import os
import argparse
from ete3 import Tree
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description = 'Plot consistency of BUSCO pairwise comparisons')
parser.add_argument('--d' , help = 'Path to BUSCO quality check pairwise comparisons')
parser.add_argument('--tree', help = 'Phylogenetic tree')
parser.add_argument('--refGenome', help = 'Name of species used as query')
parser.add_argument('--species_list', help = 'Tab delimited species list file')
parser.add_argument('--o', help = 'Output directory')




def get_species_group(species_list):
	mygroup = {}
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, phylo_class, species_name, group = line.strip().split()
			mygroup[species_name] = group
	return mygroup




def plot_consistency(mygroup, tree, refGenome, directory, output_directory):
	inconsistent_genes = defaultdict(list)
	mydict = defaultdict(list)
	t = Tree(tree)
	# Reroot the tree to the outgroup
	t.set_outgroup('Hydropsyche_tenuis')
	for node in t.traverse("postorder"):
		if node.is_leaf():
			if node.name != refGenome and node.name != "Hydropsyche_tenuis":
				num_inconsistent_genes, num_consistent_genes = number_consistent_inconsistent_genes(refGenome, node.name, directory, inconsistent_genes, mydict)

	x, y = [], []
	for key in num_consistent_genes:
		num_genes = len(num_consistent_genes[key])
		x.append(key[1])
		y.append(num_genes)

	list_colors = {'Noctuoidea': 'y', 'Bombycoidea': 'peru', 'Geometroidea': 'palevioletred', 'Drepanoidea': 'steelblue', 'Pyraloidea': 'gold', 'Papilionoidea': 'darkturquoise', 'Hesperioidea': 'darkgray',
	'Gelechioidea': 'coral', 'Zygaeinoidea': 'yellow', 'Cossoidea': 'slateblue', 'Torticoidea': 'yellowgreen', 'Tineoidea': 'cornflowerblue'}

	list_group = [mygroup[species] for species in x]
	color = [list_colors[x] for x in list_group]
	fig, ax = plt.subplots(figsize=(15, 8))
	plt.xticks(rotation = 90, ha = 'right', fontsize = 6)
	plt.bar(x, y, color = color)
	plt.ylabel('Number of consistent genes')
	plt.ylim(0, 2451)
	figure = Path(output_directory, 'number_consistent_genes.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')



	file = Path(output_directory, 'inconsistent_genes.bed')
	with open(file, 'w') as out:
		for gene in num_inconsistent_genes:
			num_species = len(num_inconsistent_genes[gene])
			out.write('{}\t{}\n'.format(gene, num_species))


def number_consistent_inconsistent_genes(query, target, directory, inconsistent_genes, mydict):
	path = Path(directory, query + '_vs_' + target)
	quality_check = list(Path(path).rglob('*.qc'))
	for file in quality_check:
		file_size = os.path.getsize(file)
		n = 0
		if file_size == 0:
			n += 1
			gene_name = Path(file).stem
			inconsistent_genes[gene_name].append(n)
		else:
			with open(file) as f:
				for line in f:
					qGenome, qChrom, qStart, qEnd, tGenome, tChrom, tStart, tEnd, number_bases = line.strip().split()
					number_bases = int(number_bases)
					mydict[query,  target].append(number_bases)
	return inconsistent_genes, mydict



if __name__ == "__main__":
	args = parser.parse_args()
	species_group = get_species_group(args.species_list)
	consistency = plot_consistency(species_group, args.tree, args.refGenome, args.d, args.o)

