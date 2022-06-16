#!/usr/bin/env python


# Author : @cb46




import argparse
import subprocess
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from collections import defaultdict




parser = argparse.ArgumentParser(description = 'Plot distribution of GERP score along the genome and in CDS/introns')
parser.add_argument('--bed', help = 'Path to species folder with GERP score bed files')
parser.add_argument('--CDS', help = 'Genomic coordinates of coding sequences')
parser.add_argument('--intron', help = 'Genomic coordinates of intron')
parser.add_argument('--o', help = 'Output directory')



def genome_wide_distribution_GERPscore(bed_f, path):
	mygerp = defaultdict(list)
	for bed in bed_f:
		species_name = '_'.join(Path(bed).stem.split('_')[0:2]).capitalize()
		with open(bed) as f:
			for line in f:
				chromosome, start, end, nucleotide, neutral_rate, rejected_substitution_score = line.strip().split()
				rejected_substitution_score = float(rejected_substitution_score)
				if rejected_substitution_score >= 0:
					if chromosome == 'W':
						mygerp[chromosome].append(rejected_substitution_score)
					elif chromosome == 'Z':
						mygerp[chromosome].append(rejected_substitution_score)
					else:
						mygerp['Autosomes'].append(rejected_substitution_score)


	color = {'Autosomes': 'red', 'W': 'gold', 'Z': 'green'}
	fig = plt.figure(figsize=(6, 4))
	# Plot density distribution of GERP score genome-wide
	for chrom in sorted(mygerp):
		scores = mygerp[chrom]
		sns.kdeplot(x=scores, fill=True, common_norm=False, alpha=.5, linewidth=0, label=chrom, color=color[chrom])
	plt.xlabel('GERP score')
	plt.ylabel('Density')
	plt.legend(prop={'size': 8})
	plt.title(species_name)
	plt.savefig(Path(path, species_name + '.genome.wide.GERPscore.pdf'), dpi = 500, bbox_inches = 'tight')
	return species_name


def distribution_GERPscore_coding_sequence_intron(cds_f, intron_f, bed_f, species_name):
	# Set empty dictionaries
	mycds, myintron = defaultdict(list), defaultdict(list)

	cds_gerp = run_bedtools(cds_f, bed_f, mycds)
	intron_gerp = run_bedtools(intron_f, bed_f, myintron)

	# Plot density distribution of GERP score in CDS and introns
	color = {'Autosomes': 'red', 'W': 'gold', 'Z': 'green'}
	fig, ax = plt.subplots(1, 2, figsize=(12, 4))
	for chrom in sorted(cds_gerp):
		scores = cds_gerp[chrom]
		sns.kdeplot(x=scores, fill=True, common_norm=False, alpha=.5, linewidth=0, label=chrom, color=color[chrom], ax = ax[0])
	ax[0].set_ylabel('Density')
	ax[0].set_xlabel('GERP score CDS')
	for chrom in sorted(intron_gerp):
		scores = intron_gerp[chrom]
		sns.kdeplot(x=scores, fill=True, common_norm=False, alpha=.5, linewidth=0, label=chrom, color=color[chrom], ax = ax[1])
	ax[1].set_xlabel('GERP score intron')
	plt.legend(prop={'size': 8})
	plt.title(species_name)
	plt.savefig(Path(path, species_name + '.genomic.regions.GERPscore.pdf'), dpi = 500, bbox_inches = 'tight')



def run_bedtools(feature_f, bed_f, dictionary):
	command = 'cat %s/*.bed | sortBed | bedtools intersect -a %s -b stdin -wa -wb' %(bed_f, feature_f)
	cmd = subprocess.check_output(command, shell = True).decode()
	outcmd = cmd.split('\n')
	for line in outcmd:
		if line:
			chrom, begin, stop, feature, strand, gene, chromosome, start, end, nucleotide, neutral_rate, rejected_substitution_score = line.strip().split()
			rejected_substitution_score = float(rejected_substitution_score)
			if rejected_substitution_score >= 0:
				if chromosome == 'W':
					dictionary[chromosome].append(rejected_substitution_score)
				elif chromosome == 'Z':
					dictionary[chromosome].append(rejected_substitution_score)
				else:
					dictionary['Autosomes'].append(rejected_substitution_score)
	return dictionary




if __name__ == "__main__":
	args = parser.parse_args()
	bed_f = list(Path(args.bed).rglob('*.bed'))
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	species_name = genome_wide_distribution_GERPscore(bed_f, path)
	distribution_cds_intron = distribution_GERPscore_coding_sequence_intron(args.CDS, args.intron, args.bed, species_name)

