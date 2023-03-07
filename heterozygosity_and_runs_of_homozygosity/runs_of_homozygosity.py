#!/usr/bin/env python



# Author : @cb46



import argparse
import pandas as pd
from pathlib import Path
from statistics import mean
from collections import defaultdict
from itertools import islice


parser = argparse.ArgumentParser(description = 'Identify stretches of homozygosity (or runs of homozygosity)')
parser.add_argument('--het', help = 'Heterozygosity file')
parser.add_argument('--w', help = 'Window size [default = 10000]', type = int, default = 10000)
parser.add_argument('--t1', help = 'Number of consecutive bins to analyse at once [default = 10]', type = int, default = 10)
parser.add_argument('--t2', help = 'Minimum number of well covered sites [default = 6000]', type = int, default = 6000)
parser.add_argument('--o', help = 'Output directory')




def filter_consecutive_bins(filename, average_genome_wide_heterozygosity, number_consecutive_bins, min_num_covered_sites):
	outsideROH = []
	insideROH = []
	lines = defaultdict(list)
	with open(filename, 'r') as infile:
		for line in infile:
			chromosome, start, end, ncov, nsites, ncorrectedHet = line.strip().split()
			start = int(start)
			end = int(end)
			ncov = int(ncov)
			ncorrectedHet = float(ncorrectedHet)
			lines[chromosome].append([start, end, ncov, ncorrectedHet])
		outsideROH, insideROH = filter_bins(lines, average_genome_wide_heterozygosity, number_consecutive_bins, min_num_covered_sites, outsideROH, insideROH)
	return outsideROH, insideROH



def filter_bins(lines, average_genome_wide_heterozygosity, number_consecutive_bins, min_num_covered_sites, outsideROH, insideROH):
	for key in lines:
		listL = lines[key]
		# Consider 10 consecutive bins at a time
		for i in range(0, len(listL), number_consecutive_bins):
			slice = list(islice(listL, i, i+number_consecutive_bins))
			# Calculate average heterozygosity within N consecutive bins
			consecutive_bins_het = [i[3] for i in slice if i[2] >= min_num_covered_sites]
			if len(consecutive_bins_het) > 1:
				het_bin = mean(consecutive_bins_het)
			else:
				het_bin = consecutive_bins_het
			if slice and consecutive_bins_het:
				for (start, end, ncov, ncorrectedHet) in slice:
					if ncov >= min_num_covered_sites:
						if het_bin <= 0.25 * average_genome_wide_heterozygosity:
							insideROH.append([key, start, end, ncov, ncorrectedHet])
						else:
							outsideROH.append([key, start, end, ncov, ncorrectedHet])
	return outsideROH, insideROH



def identify_runs_of_homozygosity(insideROH, average_genome_wide_heterozygosity, window, path, output_file):
	prev_bin = 0
	consecutive_bins = []
	sublist = []
	for (chromosome, start, end, ncov, ncorrectedHet) in insideROH:
		if ncorrectedHet <= 2 * average_genome_wide_heterozygosity:
			cur_bin = end
			if cur_bin == prev_bin + window or not sublist:
				sublist.append([chromosome, start, end, ncorrectedHet])
				prev_bin = cur_bin
			else:
				consecutive_bins.append(sublist)
				sublist = []
				sublist.append([chromosome, start, end, ncorrectedHet])
				prev_bin = cur_bin
	with open(Path(path, output_file), 'w') as output:
		for item in consecutive_bins:
			avg_het = mean([i[3] for i in item])
			if avg_het <= 0.25 * average_genome_wide_heterozygosity:
				output.write('{}\t{}\t{}\t{}\t{}\n'.format(item[0][0], item[0][1], item[-1][2], len(item),avg_het))



def save_to_output_file(outsideROH, path, output_file):
	with open(Path(path, output_file), 'w') as output:
		for item in outsideROH:
			items = '\t'.join(map(str, item))
			output.write('{}\n'.format(items))



if __name__ == "__main__":
	args = parser.parse_args()
	# Generate folder if it doesn't exist 
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	# Select bins with at least 6000 well-covered sites
	df = pd.read_csv(args.het, sep="\t", header=None, names=['Chrom', 'Start', 'End', 'Ncov', 'nHet', 'SNPcount'])
	df_f = df[df['Ncov'] >= args.t2]
	# Calculate average heterozygosity for those bins that passed the first filtering criteria
	average_genome_wide_heterozygosity = df_f['SNPcount'].mean()
	# Set output files
	output_file = Path(args.het).stem.replace('.heterozygosity', '.outsideROH.txt')
	output_file_ROH = Path(args.het).stem.replace('.heterozygosity', '.insideROH.txt')
	non_homozygous_stretches, candidate_ROH = filter_consecutive_bins(args.het, average_genome_wide_heterozygosity, args.t1, args.t2)
	runs_of_homozygosity = identify_runs_of_homozygosity(candidate_ROH, average_genome_wide_heterozygosity, args.w, path, output_file_ROH)
	save_to_output = save_to_output_file(non_homozygous_stretches, path, output_file)

  
  
