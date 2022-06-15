#!/usr/bin/env python


# Author : @cb46


import argparse
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(description = 'Plot ancestral conserved elements')
parser.add_argument('--a', help = 'Path to ancestral conserved elements files')
parser.add_argument('--g', help = 'Path to tab delimited genome files structured as follow: <chromName><TAB><chromSize>')
parser.add_argument('--species_list', help = 'Tab delimited species list file')
parser.add_argument('--l', help = 'Minimum length to retain an ancestral conserved element', type=int, default=50)
parser.add_argument('--i', help = 'Minimum identity score to retain an ancestral conserved element', type=int, default=70)
parser.add_argument('--o', help = 'Output directory')




def species_superfamily(species_list):
	mydict = {}
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, phylo_class, species_name, superfamily = line.strip().split()
			mydict[species_name] = superfamily
	return mydict



def filter_ancestral_conserved_elements(ancestral_f, min_length, min_identity):
	myelem = defaultdict(list)
	print (f"We are going to retain ancestral conserved elements that are at least {min_length} long and that have an identity score of at least {min_identity}")
	for file in ancestral_f:
		genome_name = Path(file).stem.replace('ancestral_conserved_elements_', '')
		species_name = '_'.join(genome_name.split('_')[0:2]).capitalize()
		with open(file) as f:
			for line in f:
				tgtChrom, tgtStart, tgtEnd, tgtStrand, tgtRep, tgtScore, tgtSeq, srcChrom, srcStart, srcEnd, srcLen, srcStrand, srcRep, srcSeq = line.strip().split()
				# We will only consider ancestral conserved elements outside repeats
				if float(tgtRep) == 0.0:
					# Let's further filter ancestral conserved elements based on minimum length and minimum identity score
					if int(srcLen) >= min_length and float(tgtScore) >= min_identity:
						length = int(srcLen)
						myelem[species_name].append(length)
	return myelem



def calculate_total_genome_size(genome_size_f):
	mygenome = {}
	for file in genome_size_f:
		genome_name = Path(file).stem
		species_name = '_'.join(genome_name.split('_')[0:2]).capitalize()
		g_len = 0
		with open(file) as f:
			for line in f:
				chromosome, length = line.strip().split()
				g_len += int(length)
		mygenome[species_name] = g_len
	return mygenome



def plot_ancestral_conserved_elements(myelem, mydict, mygenome, path):
	list_colors = {'Noctuoidea': 'y', 'Bombycoidea': 'peru', 'Geometroidea': 'palevioletred', 'Drepanoidea': 'steelblue', 'Pyraloidea': 'gold', 'Papilionoidea': 'darkturquoise', 'Hesperioidea': 'darkgray',
	'Gelechioidea': 'coral', 'Zygaeinoidea': 'yellow', 'Cossoidea': 'slateblue', 'Torticoidea': 'yellowgreen', 'Tineoidea': 'cornflowerblue'}
	x, y, z, colors = [], [], [], []
	for species_name in myelem:
		# Total length of ancestral conserved elements
		len_aCEs = myelem[species_name]
		total_len_aCEs = sum(len_aCEs)
		aCEs_mb = total_len_aCEs / 1000000
		# Genome size in Mb
		genome_mb = mygenome[species_name] / 1000000
		# Total number of ancestral conserved elements
		num_aCEs = len(myelem[species_name])
		# Get color for each superfamily
		superfamily = mydict[species_name]
		superfamily_color = list_colors[superfamily]
		# Append values for plot
		x.append(genome_mb)
		y.append(aCEs_mb)
		z.append(num_aCEs)
		colors.append(superfamily_color)
		len_distribution = plot_len_distribution(species_name, len_aCEs, superfamily_color, path)
	genome_vs_aCEs = plot_genome_size_vs_aCEs(x, y, colors, path)
	number_vs_length = plot_number_vs_length_aCEs(y,z, colors, path)


def plot_len_distribution(label, value, color, path):
	figure = Path(path, label+'_len_distribution.pdf')
	fig = plt.figure(figsize=(8, 6))
	plt.hist(value, density = False, bins = 50, color = color)
	plt.title(label, fontsize=8)
	plt.ylabel('Count')
	plt.xlabel('Length ancestral conserved elements (bp)')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')
	plt.clf()


def plot_genome_size_vs_aCEs(x, y, colors, path):
	figure = Path(path, 'correlation_genome_size_aCEs.pdf')
	plt.scatter(x, y, c = colors, alpha = 0.5, s = 100, edgecolor = 'black')
	plt.ylabel('Total length ancestral conserved elements (Mb)')
	plt.xlabel('Total genome length (Mb)')
	plt.ylim(0, 25)
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')
	plt.clf()


def plot_number_vs_length_aCEs(y, z, colors, path):
	figure = Path(path, 'correlation_number_vs_size_aCEs.pdf')
	plt.scatter(y, z, c = colors, alpha = 0.5, s = 100, edgecolor = 'black')
	plt.xlabel('Total length of ancestral conserved elements (Mb)')
	plt.ylabel('Total number of ancestral conserved elements')
	#plt.xlim(0, 25)
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')
	plt.clf()






if __name__ == "__main__":
	args = parser.parse_args()
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	sp_superfamily = species_superfamily(args.species_list)
	# Ancestral conserved elements
	list_a = list(Path(args.a).rglob('*.txt'))
	filtered_aCEs = filter_ancestral_conserved_elements(list_a, args.l, args.i)
	# Genome size files
	list_g = list(Path(args.g).rglob('*.bed'))
	genome_size = calculate_total_genome_size(list_g)
	# Plot
	plot_aCEs = plot_ancestral_conserved_elements(filtered_aCEs, sp_superfamily, genome_size, args.o)

	
