#!/usr/bin/env python



# Author : @cb46




import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
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
	num_roh, len_roh = [], []
	list_species, list_genome = [], []
	list_superfamily = []

	for key in sorted(myroh):
		list_species.append(key[0])
		list_superfamily.append(key[1])
		list_genome.append(key[2])
		total_roh_len = sum(myroh[key])
		total_roh_num = len(myroh[key])
		num_roh.append(total_roh_num)
		len_roh.append(total_roh_len)

		# Classify ROHs
		short_ROH = list(filter(lambda length: length <= 0.1, myroh[key]))
		medium_ROH = list(filter(lambda length: length > 0.1 and length < 1.0, myroh[key]))
		long_ROH = list(filter(lambda length: length >= 1.0, myroh[key]))
		num_short.append(len(short_ROH)), num_medium.append(len(medium_ROH)), num_long.append(len(long_ROH))
		sum_short.append(sum(short_ROH)), sum_medium.append(sum(medium_ROH)), sum_long.append(sum(long_ROH))


	list_colors = {'Noctuoidea': '#B1C968', 'Bombycoidea': '#C5A07A', 'Geometroidea': '#DB98AE', 'Drepanoidea': '#8AB1C9', 'Pyraloidea': '#ECC978', 'Papilionoidea': '#66C2A5', 'Hesperioidea': '#B3B3B3', 'Gelechioidea': '#DD927E', 
    	'Zygaeinoidea': '#FCD738', 'Cossoidea': '#BE93C6', 'Torticoidea': '#CED843', 'Tineoidea': '#979EC1'}

	df = pd.DataFrame(list(zip(list_species, list_superfamily, list_genome, num_roh, len_roh, num_short, num_medium, num_long, sum_short, sum_medium, sum_long)),
	columns = ['Species_name', 'Superfamily', 'Genome_size', 'Num_ROH', 'Length_ROH', 'Num_short', 'Num_medium', 'Num_long', 'Sum_short', 'Sum_medium', 'Sum_long'])

	fig = plt.subplots(figsize=(9, 7))
	sns.scatterplot(x = 'Length_ROH', y = 'Num_ROH', data = df, palette = list_colors, hue = 'Superfamily', linewidth=1, edgecolor = 'black', s = 140)
	plt.ylabel('Total length of ROH (Mb)')
	plt.xlabel('Total genome size (Mb)')
	plt.yticks(fontsize = 14)
	plt.xticks(fontsize = 14)
	figure = Path(path, 'fraction_genome_covered_ROH.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')

	plt.figure(figsize=(7, 14))
	sns.barplot(x = 'Num_short', y = 'Species_name', data = df, color = '#ffa600', edgecolor='black', label = 'Short: < 100 Kb')
	sns.barplot(x = 'Num_medium', y = 'Species_name', data = df, color = '#bc5090', edgecolor='black', label = 'Medium: 0.1 - 1 Mb')
	sns.barplot(x = 'Num_long', y = 'Species_name', data = df, color = '#003f5c', edgecolor='black', label = 'Long: > 1 Mb')
	plt.xlabel('Total number of ROH')
	figure = Path(path, 'number_ROH_per_size_type.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')


	plt.figure(figsize=(7, 14))
	sns.barplot(x = 'Sum_short', y = 'Species_name', data = df, color = '#ffa600', edgecolor='black', label = 'Short: < 100 Kb')
	sns.barplot(x = 'Sum_medium', y = 'Species_name', data = df, color = '#bc5090', edgecolor='black', label = 'Medium: 0.1 - 1 Mb')
	sns.barplot(x = 'Sum_long', y = 'Species_name', data = df, color = '#003f5c', edgecolor='black', label = 'Long: > 1 Mb')
	plt.xlabel('Total number of ROH')
	figure = Path(path, 'sum_ROH_per_size_type.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')





if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	species_info = get_species_genome(args.species_list)
	roh_f = list(Path(args.d).rglob('*.insideROH.txt'))
	genome_size = get_genome_size(args.genome)
	fraction_genome_in_ROH = fraction_genome_covered_by_ROH(roh_f, species_info, genome_size, path)


