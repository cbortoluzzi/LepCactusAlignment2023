#!/usr/bin/env python



# Author : @cb46


import argparse
from ete3 import Tree
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Plot genome-wide coverage')
parser.add_argument('--t', help = 'Phylogenetic tree obtained from the cactus alignment')
parser.add_argument('--l', help = 'Length of the region in base pairs [default = 1000000]', type = int, default = 1000000)
parser.add_argument('--c', help = 'Path to directory with coverage output (e.g. coverage/vanessa_atalanta_gca905147765v1)')
parser.add_argument('--f', help = 'Tab delimited species list file')
parser.add_argument('--refGenome', help = 'Name of reference genome')
parser.add_argument('--o', help = 'Output directory')



mycolors = {'Noctuoidea': '#B1C968', 'Bombycoidea': '#C5A07A', 'Geometroidea': '#DB98AE', 'Drepanoidea': '#8AB1C9', 'Pyraloidea': '#ECC978', 'Papilionoidea': '#66C2A5', 'Gelechioidea': '#DD927E', 'Zygaeinoidea': '#FCD738', 'Cossoidea': '#BE93C6', 'Torticoidea': '#CED843', 'Tineoidea': '#979EC1'}



def order_species_following_phylogenetic_tree(tree, bp, cov_f):
	cov_d = defaultdict(list)
	t = Tree(tree, format = 1)
	t.set_outgroup('tinea_trinotella_gca905220615v1')
	for node in t.traverse('postorder'):
		if node.is_leaf():
			coverageD = get_coverage(node.name, cov_f, cov_d, bp)
	return coverageD



def get_coverage(node, files, dictionary, base_pairs):
	for file in files:
		with open(file) as f:
			next(f)
			for line in f:
				query, target, lengthOfReference, percentCoverage, basesCoverage = line.strip().split()
				basesCoverage = int(basesCoverage)
				if target == node:
					coverage_alignment = round((basesCoverage/base_pairs)*100, 2)
					dictionary[target].append(coverage_alignment)
	return dictionary



def get_species_superfamily(species_list):
	superfamilies = {}
	with open(species_list) as f:
		next(f)
		for line in f:
			#tol_id, pclass, order, family, latin_name, assembly = line.strip().split()
			#superfamilies[latin_name] = family
			tol_id, pclass, order, superfamily, family, latin_name, assembly = line.strip().split()
			superfamilies[latin_name] = superfamily
	return superfamilies




def plot_coverage(coverageD, superfamilies, mycolors, refGenome, path):
	speciesL, coverageL, superfamilyL = [], [], []
	for genome in coverageD:
		species_name = '_'.join(genome.split('_')[0:2]).capitalize()
		speciesL.append(species_name)
		coverageL.append(coverageD[genome])
		superfamilyL.append(superfamilies[species_name])

	reference = '_'.join(refGenome.split('_')[0:2]).capitalize()
	fig, ax = plt.subplots(figsize=(25, 10))
	bplot = ax.boxplot(coverageL, positions=range(len(coverageL)), labels=speciesL, notch=True)
	plt.xticks(rotation = 90, ha = 'right', fontsize = 16)
	plt.yticks(fontsize = 18)
	plt.ylabel('Coverage %', fontsize = 20)
	colors = [mycolors[superfamily] for superfamily in superfamilyL]
	for artist, color in zip(bplot['boxes'], colors):
		patch = mpatches.PathPatch(artist.get_path(), color=color)
		ax.add_artist(patch)
	plt.title(reference, fontsize = 20)
	figure = Path(path, refGenome + '_coverage.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')



if __name__ == "__main__":
	args = parser.parse_args()
	# List of coverage files
	cov_f = list(Path(args.c).rglob('*.cov'))
	# Generate output directory if it doesn't extist
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	# Obtain coverage for each species, following their phylogenetic placement
	coverageD = order_species_following_phylogenetic_tree(args.t, args.l, cov_f)
	# Get the name of the superfamily each species belongs to
	superfamily = get_species_superfamily(args.f)
	# Plot coverage
	boxplot_coverage = plot_coverage(coverageD, superfamily, mycolors, args.refGenome, args.o)
	
	
