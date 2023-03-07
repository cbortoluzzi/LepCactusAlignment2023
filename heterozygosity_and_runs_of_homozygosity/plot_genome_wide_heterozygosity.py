#!/usr/bin/env python



# Author : @cb46



import argparse
from pathlib import Path
from statistics import mean, median
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Plot genome-wide heterozygosity for each Lepidoptera superfamily')
parser.add_argument('--g', help = 'Path to estimated binned heterozygosity')
parser.add_argument('--mcov', help = 'Minimum number of well-covered sites [default = 6000]', type = int, default = 6000)
parser.add_argument('--f', help = 'Tab delimited species list file')
parser.add_argument('--w', help = 'Window size used to calculate the heterozygosity [default = 10000]', type = int, default = 10000)
parser.add_argument('--o', help = 'Output directory')



mycolors = {'Noctuoidea': '#B1C968', 'Bombycoidea': '#C5A07A', 'Geometroidea': '#DB98AE', 'Drepanoidea': '#8AB1C9', 'Pyraloidea': '#ECC978', 'Papilionoidea': '#66C2A5', 'Gelechioidea': '#DD927E', 'Zygaeinoidea': '#FCD738', 'Cossoidea': '#BE93C6', 'Torticoidea': '#CED843', 'Tineoidea': '#979EC1'}


def get_species_superfamily(species_list):
	mydict = {}
	with open(species_list) as f:
		next(f)
		for line in f:
			dtol, pclass, order, superfamily, family, latin_name, assembly = line.strip().split()
			genome = latin_name.lower() + '_' +  assembly.replace('GCA_', 'gca').replace('.', 'v')
			tol = dtol.split('.')[0]
			mydict[tol] = superfamily
	return mydict


def plot_species_heterozygosity(list_files, mydict, min_cov, path, window):
	avg_sem_het = defaultdict(list)
	for file in list_files:
		tol = Path(file).name.split('.', 1)[0]
		with open(file) as f:
			for line in f:
				chrom, start, end, ncov, nhet, snp_count = line.strip().split()
				if int(ncov) >= min_cov:
					snp_count = round(float(snp_count)/window, 3)
					avg_sem_het[tol].append(snp_count)

	x, y, z= [], [], []
	for tol_id in mydict:
		if avg_sem_het[tol_id]:
			snp_count_avg = mean(avg_sem_het[tol_id])
			superfamily= mydict[tol_id]
			# Append superfamilies
			x.append(superfamily)
			# Append average heterozygosity per base pair
			y.append(snp_count_avg)
			color = mycolors[superfamily]
			# Append color based on superfamily
			z.append(color)

	df = pd.DataFrame(list(zip(x, y)), columns = ['Superfamily', 'Genome_wide_heterozygosity'])
	fig = plt.subplots(figsize=(10, 7))
	sns.boxplot(x = 'Superfamily', y = 'Genome_wide_heterozygosity', data = df, palette = mycolors)
	sns.swarmplot(x = 'Superfamily', y = 'Genome_wide_heterozygosity', data = df, color ='black', size = 10)
	plt.xticks(rotation = 90, fontsize = 16)
	plt.yticks(fontsize = 16)
	plt.xlabel('')
	plt.ylabel('Genome-wide heterozygosity', fontsize = 16)
	figure = Path(path, 'genome_wide_heterozygosity.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')




if __name__ == "__main__":
	args = parser.parse_args()
	# Generate directory if it doesn't exist
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	het_f = sorted(list(Path(args.g).rglob('*.txt')))
	list_superfamilies = get_species_superfamily(args.f)
	plot_heterozygosity = plot_species_heterozygosity(het_f, list_superfamilies, args.mcov, args.o, args.w)

	
