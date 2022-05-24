#!/usr/bin/env python



# Author : @cb46



import vcf
import argparse
from pathlib import Path
from statistics import mean
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(description = '')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--species_list', help = 'Tab delimited species file')
parser.add_argument('--o', help = 'Output directory')



def get_species_group(species_list):
	mygenome = {}
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, p_class, species_name, group = line.strip().split()
			tol = tol_id.split('.')[0]
			mygenome[tol] = [group, species_name]
	return mygenome



def plot_vcf_statistics(vcf_f, tol_vcf, mygenome, output_directory):
	vcf_reader = vcf.Reader(filename=vcf_f)
	phred_scaled_quality_score, read_depth, genotype_quality_score = [], [], []
	for record in vcf_reader:
		# let's consider only bi-allelic SNPs that pass all filtering criteria
		if not record.FILTER and record.is_snp:
			for call in record.samples:
				genotype = call['GT']
				genotype_quality = call['GQ']
				genotype_depth = call['DP']
				if genotype == '0/1':
					phred_scaled_quality_score.append(record.QUAL)
					read_depth.append(genotype_depth)
					genotype_quality_score.append(genotype_quality)
	if mygenome[tol_vcf][0] == 'Moth':
		color = '#2f6694'
	else:
		color = '#ff8000'
	fig = Path(output_directory, mygenome[tol_vcf][1] + '.PHRED.pdf')
	plot_phred_score = plot_histogram_distribution(phred_scaled_quality_score, mygenome[tol_vcf][1], 'PHRED-quality score', fig, color)
	fig = Path(output_directory, mygenome[tol_vcf][1] + '.DP.pdf')
	plot_read_depth = plot_histogram_distribution(read_depth, mygenome[tol_vcf][1], 'Read depth', fig, color)
	fig = Path(output_directory, mygenome[tol_vcf][1] + '.GQ.pdf')
	plot_genotype_quality = plot_histogram_distribution(genotype_quality_score, mygenome[tol_vcf][1], 'Genotype quality', fig, color)



def plot_histogram_distribution(values, species, xlabel, fig, color_bar):
	avg = mean(values)
	plt.hist(values, bins = 50, color = color_bar)
	plt.gca().set(title = species, ylabel = 'Frequency', xlabel = xlabel)
	plt.axvline(avg, color = 'black', linestyle = 'dashed', linewidth = 1)
	min_ylim, max_ylim = plt.ylim()
	plt.text(avg*1.1, max_ylim*0.9, 'Mean: {:.2f}'.format(avg))
	plt.savefig(fig, dpi = 500, bbox_inches = 'tight')
	plt.clf()



if __name__ == "__main__":
	args = parser.parse_args()
	sample = Path(args.vcf).name.split('.')[0]
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	species_group = get_species_group(args.species_list)
	tol_vcf = Path(args.vcf).stem.split('.')[0]
	statistics = plot_vcf_statistics(args.vcf, tol_vcf, species_group, path)

	
	
