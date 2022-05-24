#!/usr/bin/env python



# Author : @cb46



import vcf
import argparse
from pathlib import Path
from statistics import mean, stdev
import matplotlib.pyplot as plt
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Plot genome-wide heterozygosity')
parser.add_argument('--f', help = 'Path to delimited  with estimated binned-heterozygosity')
parser.add_argument('--mcov', help = 'Minimum number of well-covered sites [default = 6000]', type = int, default = 6000)
parser.add_argument('--species_list', help = 'A tab delimited species list file')
parser.add_argument('--o', help = 'Output directory')


def get_species_genome(species_list):
	mydict = {}
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, p_class, species_name, group = line.strip().split()
			tol = tol_id.split('.')[0]
			mydict[tol] = [species_name, group]
	return mydict




def plot_species_heterozygosity(list_files, mydict, min_cov, path):
	avg_sem_het = defaultdict(list)
	for file in list_files:
		species = Path(file).name.split('.', 1)[0]
		with open(file) as f:
			for line in f:
				chrom, start, end, ncov, nhet, snp_count = line.strip().split()
				if int(ncov) >= min_cov:
					snp_count = round(float(snp_count)/10000, 3)
					species_name = mydict[species][0]
					species_group = mydict[species][1]
					avg_sem_het[species_name, species_group].append(snp_count)

	x_moth, y_moth, e_moth, x_butterfly, y_butterfly, e_butterfly = [], [], [], [], [], []
	for key in sorted(avg_sem_het):
		avg = mean(avg_sem_het[key])
		stde = stdev(avg_sem_het[key])
		if key[1] == 'Moth':
			x_moth.append(key[0])
			y_moth.append(avg)
			e_moth.append(stde)
		else:
			x_butterfly.append(key[0])
			y_butterfly.append(avg)
			e_butterfly.append(stde)

	fig = plt.subplots(figsize=(12,7))
	plt.errorbar(x_moth, y_moth, yerr = e_moth, fmt='o', color='#2f6694', ecolor='#2f6694', elinewidth=1, capsize=3)
	plt.errorbar(x_butterfly, y_butterfly, yerr = e_butterfly, fmt='o', color='#ff8000', ecolor='#ff8000', elinewidth=1, capsize=3)
	plt.xticks(rotation=90, fontsize = 8)
	plt.ylabel('Number of heterozygous SNPs / bp', fontsize = 8)
	figure = Path(path, 'genome_wide_heterozygosity.pdf')
	plt.savefig(figure, dpi = 500, bbox_inches = 'tight')






if __name__ == "__main__":
	args = parser.parse_args()
	het_f = sorted(list(Path(args.f).rglob('*.txt')))
	species_name = get_species_genome(args.species_list)
	plot_heterozygosity = plot_species_heterozygosity(het_f, species_name, args.mcov, args.o)
