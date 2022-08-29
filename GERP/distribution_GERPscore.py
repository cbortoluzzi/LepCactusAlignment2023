#!/usr/bin/env python



# Author : @cb46



import argparse
import multiprocessing
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Plot distribution of GERP score along the genome, as well as in protein-coding and intronic regions')
parser.add_argument('--bed', help = 'Path to species folder with GERP score .rates file')
parser.add_argument('--CDS', help = 'Genomic coordinates of protein-coding sequences annotated with GERP score')
parser.add_argument('--intron', help = 'Genomic coordinates of intronic regions annotated with GERP score')
parser.add_argument('--o', help = 'Output directory')



def genome_wide_distribution_GERP_score(bed_f, path):
	mygerp = defaultdict(list)
	for bed in bed_f:
		species_name = '_'.join(Path(bed).stem.split('_')[0:2]).capitalize()
		with open(bed) as f:
			for line in f:
				chromosome, start, end, nucleotide, neutral_rate, rejected_substitution_score = line.strip().split()
				rejected_substitution_score = float(rejected_substitution_score)
				if chromosome == 'Z':
					mygerp[chromosome].append(rejected_substitution_score)
				elif isinstance(int(chromosome), int):
					mygerp['Autosomes'].append(rejected_substitution_score)
				elif chromosome == 'W':
					pass

	color = {'Autosomes': 'red', 'Z': 'green'}
	fig = plt.figure(figsize=(6, 4))
	# Plot density distribution of GERP score along the genome
	for chrom in sorted(mygerp):
		scores = mygerp[chrom]
		sns.kdeplot(x=scores, fill=True, common_norm=False, alpha=.5, linewidth=2, label=chrom, color=color[chrom])
	plt.xlabel('Rejected substitution score')
	plt.ylabel('Density')
	plt.legend(prop={'size': 8})
	plt.title(species_name)
	plt.savefig(Path(path, species_name + '.genome.wide.GERPscore.pdf'), dpi = 500, bbox_inches = 'tight')
	return species_name



def distribution_GERPscore_coding_sequence_intron(cds_f, intron_f, bed_f, species_name):
	# Set empty dictionaries
	mycds, myintron = defaultdict(list), defaultdict(list)

	cds_gerp = gerp_score_annotation(cds_f, mycds)
	intron_gerp = gerp_score_annotation(intron_f, myintron)

	# Plot density distribution of GERP score in CDS and introns
	color = {'Autosomes': 'red', 'Z': 'green'}
	fig, ax = plt.subplots(1, 2, figsize=(12, 4))
	for chrom in sorted(cds_gerp):
		scores = cds_gerp[chrom]
		sns.kdeplot(x=scores, fill=True, common_norm=False, alpha=.5, linewidth = 2, label=chrom, color=color[chrom], ax = ax[0])
	ax[0].set_ylabel('Density')
	ax[0].set_xlabel('GERP score CDS')
	ax[0].title.set_text(species_name)
	for chrom in sorted(intron_gerp):
		scores = intron_gerp[chrom]
		sns.kdeplot(x=scores, fill=True, common_norm=False, alpha=.5, linewidth = 2, label=chrom, color=color[chrom], ax = ax[1])
	ax[1].set_xlabel('GERP score intron')
	ax[1].title.set_text(species_name)
	plt.legend(prop={'size': 8})
	plt.savefig(Path(path, species_name + '.genomic.regions.GERPscore.pdf'), dpi = 500, bbox_inches = 'tight')



def gerp_score_annotation(file, dictionary):
	with open(file) as f:
		for line in f:
			chromosome, start, end, annotation, strand, gene, number, chrom, begin, stop, nucleotide, neutral_rate, rejected_substitution_score = line.strip().split()
			rejected_substitution_score = float(rejected_substitution_score)
			if chromosome == 'Z':
				dictionary[chromosome].append(rejected_substitution_score)
			elif isinstance(int(chromosome), int):
				dictionary['Autosomes'].append(rejected_substitution_score)
			elif chromosome == 'W':
				pass
	return dictionary



if __name__ == "__main__":
	args = parser.parse_args()
	bed_f = list(Path(args.bed).rglob('*.conserved.bed'))
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	genome_distribution_GERP = genome_wide_distribution_GERP_score(bed_f, path)
	distribution_cds_intron = distribution_GERPscore_coding_sequence_intron(args.CDS, args.intron, args.bed, genome_distribution_GERP)


