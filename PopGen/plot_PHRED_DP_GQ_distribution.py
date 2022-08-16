#!/usr/bin/env python



# Author : @cb46



import vcf
import argparse
from pathlib import Path
from statistics import mean
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(description = 'Plot distribution of phred-quality score, read depth, and genotype quality')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--species_name', help = 'Species name')
parser.add_argument('--o', help = 'Output directory')



def plot_vcf_statistics(vcf_f, species_name, output_file, output_directory):
	vcf_reader = vcf.Reader(filename=vcf_f)
	phred_scaled_quality_score, read_depth, genotype_quality_score = [], [], []
	for record in vcf_reader:
		# Let's consider only bi-allelic SNPs
		if not record.FILTER and record.is_snp:
			for call in record.samples:
				genotype = call['GT']
				genotype_quality = call['GQ']
				genotype_depth = call['DP']
				# Let's consider only heterozygous sites
				if genotype == '0/1':
					phred_scaled_quality_score.append(record.QUAL)
					read_depth.append(genotype_depth)
					genotype_quality_score.append(genotype_quality)

	# Distribution of phred-quality score
	fig = Path(output_directory, output_file + '.QUAL.pdf')
	plot_phred_score = plot_histogram_distribution(phred_scaled_quality_score, species_name, 'PHRED-quality score', fig)
	# Distribution of read depth
	fig = Path(output_directory, output_file + '.DP.pdf')
	plot_read_depth = plot_histogram_distribution(read_depth, species_name, 'Read depth', fig)
	# Distribution of genotype quality
	fig = Path(output_directory, output_file + '.GQ.pdf')
	plot_genotype_quality = plot_histogram_distribution(genotype_quality_score, species_name, 'Genotype quality', fig)



def plot_histogram_distribution(values, species, xlabel, fig):
	avg = mean(values)
	plt.hist(values, bins = 50)
	plt.gca().set(title = species, ylabel = 'Frequency', xlabel = xlabel)
	plt.axvline(avg, color = 'black', linestyle = 'dashed', linewidth = 1)
	min_ylim, max_ylim = plt.ylim()
	plt.text(avg*1.1, max_ylim*0.9, 'Mean: {:.2f}'.format(avg))
	plt.savefig(fig, dpi = 500, bbox_inches = 'tight')
	plt.clf()



if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	output_file = Path(args.vcf).stem
	statistics = plot_vcf_statistics(args.vcf, args.species_name, output_file, path)
	
	
