#!/usr/bin/env python


# Author : @cb46


import os
import argparse
from ete3 import Tree
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description = 'Plot consistency of BUSCO/OrthoFinder pairwise comparisons')
parser.add_argument('--d' , help = 'Path to BUSCO/OrthoFinder quality check pairwise comparisons')
parser.add_argument('--tree', help = 'Phylogenetic tree used as guide tree in cactus')
parser.add_argument('--refGenome', help = 'Name of species used as query')
parser.add_argument('--species_list', help = 'Tab delimited species list file')
parser.add_argument('--o', help = 'Output directory')




def get_species_group(species_list):
	superfamilies = {}
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, phylo_class, species_name, superfamily = line.strip().split()
			superfamilies[species_name] = superfamily
	return superfamilies




def plot_consistency(superfamilies, tree, refGenome, directory, output_directory, label):
	inconsistent_genes = defaultdict(list)
	mydict = {}
	t = Tree(tree)
	# Reroot the tree to the outgroup
	t.set_outgroup('Hydropsyche_tenuis')
	for node in t.traverse("postorder"):
		if node.is_leaf():
			# Exclude outgroup
			if node.name != refGenome and node.name != "Hydropsyche_tenuis":
				num_inconsistent_genes, num_consistent_genes = number_consistent_inconsistent_genes(refGenome, node.name, directory, inconsistent_genes, mydict)

	x, y = [], []
	for key in num_consistent_genes:
		num_genes = num_consistent_genes[key]
		if num_genes > 0:
			x.append(key[1])
			y.append(num_genes)

	list_colors = {'Noctuoidea': '#B1C968', 'Bombycoidea': '#C5A07A', 'Geometroidea': '#DB98AE', 'Drepanoidea': '#8AB1C9', 'Pyraloidea': '#ECC978', 'Papilionoidea': '#66C2A5', 'Hesperioidea': '#B3B3B3', 'Gelechioidea': '#DD927E', 
        'Zygaeinoidea': '#FCD738', 'Cossoidea': '#BE93C6', 'Torticoidea': '#CED843', 'Tineoidea': '#979EC1'}
	list_superfamilies = [superfamilies[species] for species in x]
	color = [list_colors[x] for x in list_superfamilies]
	fig, ax = plt.subplots(figsize=(15, 8))
	plt.xticks(rotation = 90, ha = 'right', fontsize = 12)
	plt.bar(x, y, color = color, edgecolor = 'black')
	plt.ylabel('Number of consistent ' + label)
	figure = Path(output_directory, 'number_consistent_' + label + '.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')

	file = Path(output_directory, 'inconsistent_' + label + '.bed')
	with open(file, 'w') as out:
		for gene in num_inconsistent_genes:
			num_species = len(num_inconsistent_genes[gene])
			out.write('{}\t{}\n'.format(gene, num_species))


def number_consistent_inconsistent_genes(query, target, directory, inconsistent_genes, mydict):
	path = Path(directory, query + '_vs_' + target)
	quality_check = list(Path(path).rglob('*.qc'))
	count = 0
	for file in quality_check:
		file_size = os.path.getsize(file)
		n = 0
		# A gene/orthogroup with a file size of 0 is classified as inconsistent
		if file_size == 0:
			n += 1
			gene_name = Path(file).stem
			inconsistent_genes[gene_name].append(n)
		# Otherwise as consistent
		else:
			with open(file) as f:
				for line in f:
					qGenome, qChrom, qStart, qEnd, tGenome, tChrom, tStart, tEnd, number_bases = line.strip().split()
					count += 1
	mydict[query,  target] = count
	return inconsistent_genes, mydict



if __name__ == "__main__":
	args = parser.parse_args()
	# Generate folder if it doesn't exist
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	if 'orthofinder' in args.d:
		label = 'orthogroups'
	else:
		label = 'BUSCO'
	species_superfamily = get_species_group(args.species_list)
	consistency = plot_consistency(species_superfamily, args.tree, args.refGenome, args.d, args.o, label)

	
