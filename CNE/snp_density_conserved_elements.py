#!/usr/bin/env python



# Author : @cb46



import vcf
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Estimate SNP density inside and outside lepidoptera-specific ancestral conserved elements')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--a', help = 'Ancestral conserved elements')
parser.add_argument('--o', help = 'Output directory')




def regions_outside_conserved_elements(ancestral_f):
	genome = {}
	mydict = defaultdict(list)
	final = {}
	with open(ancestral_f) as f:
		for line in f:
			chromosome, start, length, strand, chromosome_len, perc_identity, perc_repeats, sequence, element_type = line.strip().split()
			start = int(start)
			end = start + int(length)
			chromosome_len = int(chromosome_len)
			genome[chromosome] = chromosome_len
			mydict[chromosome].append([start, end])

	for chrom in mydict:
		list_starts = sorted([i[0] for i in mydict[chrom]])
		list_ends = sorted([i[1] for i in mydict[chrom]])
		if list_ends[-1] != genome[chrom]:
			start = list_ends[-1]
			end = genome[chrom]
			final[chrom, start, end] = 'not conserved'
		if list_starts[0] != 0:
			final[chrom, 0, list_starts[0]] = 'not conserved'
		regions_in_between = list(zip(list_ends, list_starts[1:]))
		for (start, end) in regions_in_between:
			final[chrom, start, end] = 'not conserved'
	return final



def count_heterozygous_sites_inside(vcf_f, ancestral_f, output_file, path):
	snp_density = {}
	with open(ancestral_f) as f:
		for line in f:
			chromosome, start, length, strand, chromosome_len, perc_identity, perc_repeats, sequence, element_type = line.strip().split()
			start = int(start)
			end = start + int(length)
			try:
				snp_density = estimate_snp_density(chromosome, start, end, vcf_f, snp_density)
			except ValueError:
				continue
	with open(Path(path, output_file), 'w') as out:
		for key, value in snp_density.items():
			out.write('{}\t{}\t{}\t{}\t{}\n'.format(key[0], key[1], key[2], value, 'conserved_element'))



def count_heterozygous_sites_outside(vcf_f, final, output_file, path):
	snp_density = {}
	for key in sorted(final.keys()):
		start = key[1]
		end = key[2]
		try:
			snp_density = estimate_snp_density(key[0], start, end, vcf_f, snp_density)
		except ValueError:
			continue
	with open(Path(path, output_file), 'a') as out:
		for i, j in snp_density.items():
			out.write('{}\t{}\t{}\t{}\t{}\n'.format(i[0], i[1], i[2], j, 'not_conserved'))



def estimate_snp_density(chromosome, start, end, vcf_f, dictionary):
	vcf_reader = vcf.Reader(filename=vcf_f)
	nhet = 0
	for record in vcf_reader.fetch(chromosome, start, end):
		nhet += record.num_het
	dictionary[chromosome, start, end] = nhet
	return dictionary





if __name__ == "__main__":
	args = parser.parse_args()
	# Create output directory if it doesn't exist
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	output_file = Path(args.a).stem + '.SNPdensity.txt'
	conserved_elements_outside = regions_outside_conserved_elements(args.a)
	snp_density_inside = count_heterozygous_sites_inside(args.vcf, args.a, output_file, path)
	snp_density_outside = count_heterozygous_sites_outside(args.vcf, conserved_elements_outside, output_file, path)

  
