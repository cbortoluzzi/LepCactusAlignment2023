#!/usr/bin/env python



# Author : @cb46



import vcf
import argparse
from pathlib import Path
from statistics import mean
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(description = '')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--species' , help = 'Species name')
parser.add_argument('--group' , help = 'Group the species belongs to')
parser.add_argument('--o', help = 'Output directory')



def plot_vcf_statistics(vcf_f, sample, output_directory, species, color):
	vcf_reader = vcf.Reader(filename=vcf_f)
	phred_scaled_quality_score, read_depth, genotype_quality_score = [], [], []
	for record in vcf_reader:
		# let's consider only bi-allelic SNPs that pass all filtering criteria
		if not record.FILTER and record.is_snp:
			for call in record.samples:
				genotype = call['GT']
				genotype_quality = call['GQ']
				genotype_depth = call['DP']
				if genotype == '0/1' or genotype == '1/1':
					phred_scaled_quality_score.append(record.QUAL)
					read_depth.append(genotype_depth)
					genotype_quality_score.append(genotype_quality)

	fig = Path(output_directory, sample + '.qual.pdf')
	plot_phred_score = plot_histogram_distribution(phred_scaled_quality_score, species, 'PHRED-quality score', fig, color)
	fig = Path(output_directory, sample + '.dp.pdf')
	plot_read_depth = plot_histogram_distribution(read_depth, species, 'Read depth', fig, color)
	fig = Path(output_directory, sample + '.gq.pdf')
	plot_genotype_quality = plot_histogram_distribution(genotype_quality_score, species, 'Genotype quality', fig, color)



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
	if args.group == 'Moth' or args.group == 'moth':
		color = '#2f6694'
	else:
		color = '#ff8000'
	statistics = plot_vcf_statistics(args.vcf, sample, path, args.species, color)
