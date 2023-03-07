#!/usr/bin/env python



# Author : @cb46




import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Plot runs of homozygosity (ROHs) along the genome')
parser.add_argument('--genome', help = 'Tab delimited file with genome size information for each species, one per line')
parser.add_argument('--d', help = 'Path to folder with runs of homozygosity')
parser.add_argument('--f', help = 'A tab delimited species list file')
parser.add_argument('--o', help = 'Output directory')



mycolors = {'Noctuoidea': '#B1C968', 'Bombycoidea': '#C5A07A', 'Geometroidea': '#DB98AE', 'Drepanoidea': '#8AB1C9', 'Pyraloidea': '#ECC978', 'Papilionoidea': '#66C2A5', 'Gelechioidea': '#DD927E', 'Zygaeinoidea': '#FCD738', 'Cossoidea': '#BE93C6', 'Torticoidea': '#CED843', 'Tineoidea': '#979EC1'}



def get_species_information(species_list):
	mydict = {}
	with open(species_list) as f:
		for line in f:
			tol_id, pclass, order, superfamily, family, latin_name, assembly  = line.strip().split()
			tol = tol_id.split('.')[0]
			mydict[tol] = [latin_name, superfamily]
	return mydict



def get_genome_size(genome_f):
	mygen = {}
	with open(genome_f) as f:
		for line in f:
			species_name, genome_size = line.strip().split()
			# Obtain the genome size in Mb
			genome_size = int(genome_size)/1000000
			mygen[species_name] = genome_size
	return mygen



def fraction_genome_covered_by_ROH(roh_files, mydict, mygen, path):
	myroh = defaultdict(list)
	for roh_file in roh_files:
		tol = Path(roh_file).stem.split('.')[0]
		species_name = mydict[tol][0]
		superfamily = mydict[tol][1]
		genome = mygen[species_name]
		with open(roh_file) as f:
			for line in f:
				chromosome, start, end, size, avg_snp_count = line.strip().split()
				length = int(end) - int(start)
				length = length / 1000000
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

		# Based on their length, classify ROHs into short, medium, and long
		short_ROH = list(filter(lambda length: length <= 0.1, myroh[key]))
		medium_ROH = list(filter(lambda length: length > 0.1 and length < 1.0, myroh[key]))
		long_ROH = list(filter(lambda length: length >= 1.0, myroh[key]))
		num_short.append(len(short_ROH)), num_medium.append(len(medium_ROH)), num_long.append(len(long_ROH))
		sum_short.append(sum(short_ROH)), sum_medium.append(sum(medium_ROH)), sum_long.append(sum(long_ROH))

	df_num = pd.DataFrame(list(zip(list_species, list_superfamily, list_genome, num_roh, len_roh, num_short, num_medium, num_long)),
	columns = ['Species_name', 'Superfamily', 'Genome_size', 'Num_ROH', 'Length_ROH', 'Num_short', 'Num_medium', 'Num_long'])

	df_sum = pd.DataFrame(list(zip(list_species, list_superfamily, list_genome, num_roh, len_roh, sum_short, sum_medium, sum_long)),
	columns = ['Species_name', 'Superfamily', 'Genome_size', 'Num_ROH', 'Length_ROH', 'Sum_short', 'Sum_medium', 'Sum_long'])

	# Plot correlation between length and number of ROHs
	fig = plt.subplots(figsize=(9, 7))
	sns.scatterplot(x = 'Length_ROH', y = 'Num_ROH', data = df_num, palette = mycolors, hue = 'Superfamily', linewidth=1, edgecolor = 'black', s = 140)
	plt.ylabel('Total length of ROH (Mb)')
	plt.xlabel('Total genome size (Mb)')
	plt.yticks(fontsize = 14)
	plt.xticks(fontsize = 14)
	figure = Path(path, 'fraction_genome_covered_ROH.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')

	# Plot the number of ROHs as a stacked barplot
	plt.figure(figsize=(20, 7))
	plt.bar(df_num['Species_name'], df_num['Num_short'], color = '#ffa600', edgecolor='black', label = 'Short: < 100 Kb', bottom = df_num['Num_medium'])
	plt.bar(df_num['Species_name'], df_num['Num_medium'], color = '#bc5090', edgecolor='black', label = 'Medium: 0.1 - 1 Mb', bottom = df_num['Num_long'])
	plt.bar(df_num['Species_name'], df_num['Num_long'], color = '#003f5c', edgecolor='black', label = 'Long: > 1 Mb')
	plt.ylim(0, 1600)
	plt.xticks(fontsize = 14, rotation = 90)
	plt.yticks(fontsize = 14)
	plt.ylabel('Total number of ROH', fontsize = 18)
	plt.legend()
	figure = Path(path, 'number_ROH_per_size_type.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')

	# Plot the length of ROHs as a stacked barplot
	plt.figure(figsize=(20, 7))
	plt.bar(df_sum['Species_name'], df_sum['Sum_short'], color = '#ffa600', edgecolor='black', label = 'Short: < 100 Kb', bottom = df_sum['Sum_medium'])
	plt.bar(df_sum['Species_name'], df_sum['Sum_medium'], color = '#bc5090', edgecolor='black', label = 'Medium: 0.1 - 1 Mb', bottom = df_sum['Sum_long'])
	plt.bar(df_sum['Species_name'], df_sum['Sum_long'], color = '#003f5c', edgecolor='black', label = 'Long: > 1 Mb')
	plt.xticks(fontsize = 14, rotation = 90)
	plt.yticks(fontsize = 14)
	plt.ylabel('Total number of ROH', fontsize = 18)
	plt.legend()
	figure = Path(path, 'sum_ROH_per_size_type.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')




if __name__ == "__main__":
	args = parser.parse_args()
	# Generate directory if it doesn't exist
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	# List of ROHs
	roh_files = list(Path(args.d).rglob('*.insideROH.txt'))
	species_info = get_species_information(args.f)
	genome_size = get_genome_size(args.genome)
	fraction_genome_in_ROH = fraction_genome_covered_by_ROH(roh_files, species_info, genome_size, path)

  
