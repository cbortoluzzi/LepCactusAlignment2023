#!/usr/bin/env python


# Author : @cb46


import os
import argparse
from ete3 import Tree
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description = 'Plot consistency of BUSCO/OrthoFinder pairwise comparisons')
parser.add_argument('--d' , help = 'Path to BUSCO quality check pairwise comparisons')
parser.add_argument('--tree', help = 'Phylogenetic tree used as guide tree in cactus')
parser.add_argument('--refGenome', help = 'Name of species used as query')
parser.add_argument('--o', help = 'Output directory')



def phylogenetic_distance(tree, refGenome):
	distance = {}
	t = Tree(tree)
	# Reroot the tree to the outgroup
	t.set_outgroup('Hydropsyche_tenuis')
	for node in t.traverse("postorder"):
		if node.is_leaf():
			if node.name != refGenome and node.name != "Hydropsyche_tenuis":
				d = t.get_distance(refGenome, node.name)
				distance[node.name] = float(d)
	return distance



def plot_coverage_number_consistent_genes(directory, distance, path):
	mydict = defaultdict(list)
	inconsistent_genes = defaultdict(list)
	for subdir, dirs, files in os.walk(directory):
		for file in files:
			query, target = Path(subdir).stem.split('_vs_')
			qc_file = os.path.join(subdir, file)
			if qc_file.endswith('.qc'):
				file_size = os.path.getsize(qc_file)
				n = 0
				if file_size == 0:
					n += 1
					gene_name = Path(qc_file).stem
					inconsistent_genes[gene_name].append(n)
				else:
					with open(qc_file) as f:
						for line in f:
							qGenome, qChrom, qStart, qEnd, tGenome, tChrom, tStart, tEnd, number_bases = line.strip().split()
							number_bases = int(number_bases)
							p_distance = distance[target]
							mydict[p_distance, target].append(number_bases)

	x, y = [], []
	for key in sorted(mydict):
		print (key)
		num_consistent_genes = len(mydict[key])
		x.append(key[1])
		y.append(num_consistent_genes)

	fig, ax = plt.subplots(figsize=(15, 8))
	plt.xticks(rotation = 90, ha = 'right', fontsize = 6)
	plt.bar(x, y, color = '#2ca02c')
	plt.ylabel('Number of consistent genes')
	plt.ylim(0, 2500)
	plt.axhline(y=2451, color='r', linestyle='-')
	figure = Path(path, 'number_consistent_genes.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')
	plt.clf()


	file = Path(path,' inconsistent_genes.bed')
	with open(file, 'w') as out:
		for gene in inconsistent_genes:
			num_species = len(inconsistent_genes[gene])
			out.write('{}\t{}\n'.format(gene, num_species))



if __name__ == "__main__":
	args = parser.parse_args()
	p_distance = phylogenetic_distance(args.tree, args.refGenome)
	plot_genes = plot_coverage_number_consistent_genes(args.d, p_distance, args.o)



