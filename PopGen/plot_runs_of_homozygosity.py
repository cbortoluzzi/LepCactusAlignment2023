#!/usr/bin/env python



# Author : @cb46



import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Plot fraction genome covered by ROOHs')
parser.add_argument('--d', help = 'Path to folder with runs of homozygosity')
parser.add_argument('--genome', help = 'Tab delimited file with genome size information for each species, one per line')
parser.add_argument('--species_list', help = 'A tab delimited species list file')
parser.add_argument('--o', help = 'Output directory')



def get_species_genome(species_list):
	mydict = {}
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, p_class, species_name, superfamily  = line.strip().split()
			tol = tol_id.split('.')[0]
			mydict[tol] = [species_name, superfamily]
	return mydict



def get_genome_size(genome_f):
	mygen = {}
	with open(genome_f) as f:
		for line in f:
			species_name, genome_size = line.strip().split()
			genome_size = int(genome_size)/1000000
			mygen[species_name] = genome_size
	return mygen



def fraction_genome_covered_by_ROH(roh_f, mydict, mygen, path):
	myroh = defaultdict(list)
	for roh in roh_f:
		tol = Path(roh).stem.split('.')[0]
		species_name = mydict[tol][0]
		superfamily = mydict[tol][1]
		genome = mygen[species_name]
		with open(roh) as f:
			for line in f:
				chromosome, start, end, length, avg_snp_count = line.strip().split()
				length = int(length)/1000000
				myroh[species_name, superfamily, genome].append(length)

	num_short, num_medium, num_long = [], [], []
	sum_short, sum_medium, sum_long = [], [], []
	list_species, num_roh, len_roh = [], [], []
	list_superfamily, list_genome = [], []
	for key in sorted(myroh):
		list_species.append(key[0])
		total_roh_len = sum(myroh[key])
		total_roh_num = len(myroh[key])
		num_roh.append(total_roh_num)
		len_roh.append(total_roh_len)
		list_superfamily.append(key[1])
		list_genome.append(key[2])

		# Classify ROHs
		short_ROH = list(filter(lambda length: length <= 0.1, myroh[key]))
		medium_ROH = list(filter(lambda length: length > 0.1 and length < 1.0, myroh[key]))
		long_ROH = list(filter(lambda length: length >= 1.0, myroh[key]))
		num_short.append(len(short_ROH)), num_medium.append(len(medium_ROH)), num_long.append(len(long_ROH))
		sum_short.append(sum(short_ROH)), sum_medium.append(sum(medium_ROH)), sum_long.append(sum(long_ROH))

	list_colors = {'Noctuoidea': 'y', 'Bombycoidea': 'peru', 'Geometroidea': 'palevioletred', 'Drepanoidea': 'steelblue', 'Pyraloidea': 'gold', 'Papilionoidea': 'darkturquoise', 'Hesperioidea': 'darkgray',
	'Gelechioidea': 'coral', 'Zygaeinoidea': 'yellow', 'Cossoidea': 'slateblue', 'Torticoidea': 'yellowgreen', 'Tineoidea': 'cornflowerblue'}

	colors = [list_colors[x] for x in list_superfamily]
	plt.tick_params(labelsize = 8)
	plt.scatter(list_genome, len_roh, c = colors, alpha=0.5, s=100, edgecolor = 'black')
	plt.ylabel('Total length of ROH (Mb)')
	plt.xlabel('Total genome size (Mb)')
	figure = Path(path, 'fraction_genome_covered_ROH.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')
	plt.clf()

	plt.scatter(len_roh, num_roh, c = colors, alpha=0.5, s=100, edgecolor = 'black')
	plt.ylabel('Total number of ROH')
	plt.xlabel('Total length of ROH (Mb)')
	figure = Path(path, 'total_length_vs_total_number_ROH.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')
	plt.clf()

	plt.figure(figsize=(6, 8))
	plt.tick_params(labelsize = 6)
	w = 1.0
	plt.barh(list_species, num_short, w, color = '#ffa600', edgecolor='black', label = 'Short: < 100 Kb')
	plt.barh(list_species, num_medium, w, color = '#bc5090', edgecolor='black', label = 'Medium: 0.1 - 1 Mb')
	plt.barh(list_species, num_long, w, color = '#003f5c', edgecolor='black', label = 'Long: > 1 Mb')
	plt.xlabel('Total number of ROH')
	plt.legend(fontsize = 8)
	figure = Path(path, 'number_ROH_per_size_type.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')
	plt.clf()

	plt.figure(figsize=(6, 8))
	plt.tick_params(labelsize = 6)
	w = 1.0
	plt.barh(list_species, sum_short, w, color = '#ffa600', edgecolor='black', label = 'Short: < 100 Kb')
	plt.barh(list_species, sum_medium, w, color = '#bc5090', edgecolor='black', label = 'Medium: 0.1 - 1 Mb')
	plt.barh(list_species, sum_long, w, color = '#003f5c', edgecolor='black', label = 'Long: > 1 Mb')
	plt.xlabel('Total length of ROH (Mb)')
	plt.legend(fontsize = 8)
	figure = Path(path, 'sum_ROH_per_size_type.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')
	plt.clf()




if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	species_info = get_species_genome(args.species_list)
	roh_f = list(Path(args.d).rglob('*.insideROH.txt'))
	genome_size = get_genome_size(args.genome)
	fraction_genome_in_ROH = fraction_genome_covered_by_ROH(roh_f, species_info, genome_size, path)

  
  
