#!/usr/bin/env python



# Author : @cb46



import vcf
import argparse
from ete3 import Tree
from pathlib import Path
from statistics import mean, stdev
import matplotlib.pyplot as plt
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Plot genome-wide heterozygosity')
parser.add_argument('--f', help = 'Path to estimated binned-heterozygosity files')
parser.add_argument('--tree', help = 'Phylogenetic tree')
parser.add_argument('--mcov', help = 'Minimum number of well-covered sites [default = 6000]', type = int, default = 6000)
parser.add_argument('--species_list', help = 'A tab delimited species list file')
parser.add_argument('--o', help = 'Output directory')



def order_species_by_phylo(tree, species_list):
	phylo = {}
	t = Tree(tree)
	# Reroot the tree to the outgroup
	t.set_outgroup('Hydropsyche_tenuis')
	for node in t.traverse('postorder'):
		if node.is_leaf() and node.name != "Hydropsyche_tenuis":
			tree_d = get_species_genome(species_list, node.name, phylo)
	return tree_d


def get_species_genome(species_list, node, phylo):
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, p_class, species_name, superfamily = line.strip().split()
			tol = tol_id.split('.')[0]
			if species_name == node:
				phylo[tol] = [node, superfamily]
	return phylo


def plot_species_heterozygosity(list_files, min_cov, path, tree_d):
	avg_sem_het = defaultdict(list)
	for file in list_files:
		tol = Path(file).name.split('.', 1)[0]
		with open(file) as f:
			for line in f:
				chrom, start, end, ncov, nhet, snp_count = line.strip().split()
				if int(ncov) >= min_cov:
					snp_count = round(float(snp_count)/10000, 3)
					avg_sem_het[tol].append(snp_count)


	list_colors = {'Noctuoidea': '#B1C968', 'Bombycoidea': '#C5A07A', 'Geometroidea': '#DB98AE', 'Drepanoidea': '#8AB1C9', 'Pyraloidea': '#ECC978', 'Papilionoidea': '#66C2A5', 'Hesperioidea': '#B3B3B3', 'Gelechioidea': '#DD927E', 
        'Zygaeinoidea': '#FCD738', 'Cossoidea': '#BE93C6', 'Torticoidea': '#CED843', 'Tineoidea': '#979EC1'}
	
	x, y, yerr, color = [], [], [], []
	for tol_id in tree_d:
		if tol_id in avg_sem_het.keys():
			snp_count_stdev = stdev(avg_sem_het[tol_id])
			snp_count_avg = mean(avg_sem_het[tol_id])
			species_name = tree_d[tol_id][0]
			superfamily = tree_d[tol_id][1]
			x.append(species_name)
			y.append(snp_count_avg)
			yerr.append(snp_count_stdev)
			color.append(superfamily)

	colors = [list_colors[x] for x in color]
	fig = plt.subplots(figsize=(12, 7))
	plt.errorbar(x, y, yerr = yerr, fmt='o', color='black', ecolor=colors, elinewidth=2, capsize=0)
	plt.xticks(rotation=90, fontsize = 8)
	plt.ylabel('Number of corrected heterozygous SNPs / bp')
	figure = Path(path, 'genome_wide_heterozygosity.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')





if __name__ == "__main__":
	args = parser.parse_args()
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	het_f = sorted(list(Path(args.f).rglob('*.txt')))
	phylo = order_species_by_phylo(args.tree, args.species_list)
	plot_heterozygosity = plot_species_heterozygosity(het_f, args.mcov, args.o, phylo)



