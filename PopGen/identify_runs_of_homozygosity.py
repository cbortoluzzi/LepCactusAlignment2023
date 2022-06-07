#!/usr/bin/env python



# Author : @cb46



import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from statistics import mean


parser = argparse.ArgumentParser(description = 'Identify runs of homozygosity')
parser.add_argument('--het', help = 'Heterozygosity file')
parser.add_argument('--t1', help = 'Minimum number of well-covered sites', type = int)
parser.add_argument("--t2", help = "Threshold to filter the snp count in a bin [default = 0.25]", type=float, default = 0.25)
parser.add_argument('--o', help = 'Output directory')




def select_consecutive_bins(het_f, genome_wide_heterozygosity, threshold1, threshold2, path, output_f):
	prev_bin = 0
	consecutive_bins = []
	sublist = []
	outside_ROH = Path(path, output_f + '.outsideROH.txt')
	inside_ROH = Path(path, output_f + '.insideROH.txt')
	with open(het_f) as f, open(outside_ROH, 'w') as outside:
		for line in f:
			chromosome, start, end, ncov, nhet, snpcount = line.strip().split()
			ncov, snpcount = int(ncov), float(snpcount)
			if ncov >= threshold1:
				if snpcount < threshold2 * genome_wide_heterozygosity:
					cur_bin = int(end)
					if cur_bin == prev_bin + 10000 or not sublist:
						sublist.append([chromosome, start, end, snpcount])
						prev_bin = cur_bin
					else:
						consecutive_bins.append(sublist)
						sublist = []
						sublist.append([chromosome, start, end, snpcount])
						prev_bin = cur_bin
				else:
					outside.write('{}\t{}\t{}\t{}\t{}\n'.format(chromosome, start, end, ncov, snpcount))
	with open(inside_ROH, 'w') as inside:
		for item in consecutive_bins:
			avg_snp_count = round(mean([i[-1] for i in item]), 2)
			roh_start, roh_end = int(item[0][1]), int(item[-1][2])
			len_roh = roh_end - roh_start
			inside.write('{}\t{}\t{}\t{}\t{}\n'.format(item[0][0], roh_start, roh_end, len_roh, avg_snp_count))




if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	output_f = Path(args.het).stem.replace('.heterozygosity', '')
	df = pd.read_csv(args.het, sep="\t", header=None, names=['Chrom', 'Start', 'End', 'Ncov', 'nHet', 'SNPcount'])
	df_f = df[df['Ncov'] >= args.t1]
	genome_wide_heterozygosity = df_f['SNPcount'].mean()
	consecutive_bins = select_consecutive_bins(args.het, genome_wide_heterozygosity, args.t1, args.t2, path, output_f)

