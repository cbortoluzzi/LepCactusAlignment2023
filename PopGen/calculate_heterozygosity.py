#!/usr/bin/env python



# Author : @cb46



import os
import vcf
import argparse
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(description = 'Calculate genome-wide heterozygosity using a sliding window approach')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--bam', help = 'BAM file')
parser.add_argument('--cov', help = 'Genome-wide coverage')
parser.add_argument('--w', help = 'Window size [default = 10000 bp]', type = int, default = 10000)
parser.add_argument('--q', help = 'Minimum phred-quality score [default = 15]', type = int, default = 15)
parser.add_argument('--dp', help = 'Minimum read depth to retain a variant [default = 6]', type = int, default = 6)
parser.add_argument('--gq', help = 'Minimum genotype quality to retain a variant [default = 20]' , type = int, default = 20)
parser.add_argument('--o', help = 'Output directory')




def sequences_bam(bam_f):
	mygenome = {}
	# Get BAM index stats
	command = 'samtools idxstats %s | cut -f 1,2' %(bam_f)
	cmd = subprocess.check_output(command, shell = True).decode()
	outcmd = cmd.split('\n')
	for line in outcmd:
		if line:
			chromosome, length = line.strip().split()
			length = int(length)
			try:
				if str(chromosome) == 'W' or str(chromosome) == 'Z' or isinstance(int(chromosome), int):
					mygenome[chromosome] = length
			except ValueError:
				# We exclude contigs/scaffolds
				print (f"Contig {chromosome} is not considered")
	return mygenome



def average_genome_coverage(coverage):
	with open(coverage) as f:
		for line in f:
			line = line.split('=')
			avg_genome_coverage = line[1].replace(' ','')
			max_depth = 2 * float(avg_genome_coverage)
	return max_depth



def calculate_binned_heterozygosity(mygenome, window, min_depth, max_depth, bam_f, vcf_f, min_gq, min_qual, filename, path):
	for chromosome in mygenome:
		seq_length = mygenome[chromosome]
		for i in range(0, seq_length, window):
			start = i
			end = i + window
			bam_depth(chromosome, start, end, min_depth, max_depth, bam_f, vcf_f, min_gq, min_qual, window, filename, path)



def bam_depth(chromosome, start, end, min_depth, max_depth, bam_f, vcf_f, min_gq, min_qual, window, filename, path):
	cov_sites = 0
	command = 'samtools depth -r %s:%d-%d %s' %(chromosome, start, end, bam_f)
	cmd = subprocess.check_output(command, shell = True).decode()
	outcmd = cmd.split('\n')
	for line in outcmd:
		if line:
			chromosome, position, depth = line.strip().split()
			if int(depth) >= min_depth and int(depth) <= max_depth:
				cov_sites += 1
	heterozygosity(chromosome, start, end, cov_sites, vcf_f, min_depth, max_depth, min_gq, min_qual, window, filename, path)



def heterozygosity(chromosome, start, end, cov_sites, vcf_f, min_depth, max_depth, min_gq, min_qual, window, filename, path):
	vcf_reader = vcf.Reader(filename=vcf_f)
	nhet = 0
	for record in vcf_reader.fetch(chromosome, start, end):
		# let's consider only bi-allelic SNPs that pass all filtering criteria
		if not record.FILTER and record.is_snp:
			if record.QUAL >= min_qual:
				for call in record.samples:
					read_depth = call['DP']
					genotype_quality = call['GQ']
					if read_depth >= min_depth and read_depth <= max_depth and genotype_quality >= min_gq:
						nhet += record.num_het
	if cov_sites == window + 1:
		cov_sites = window
	try:
		SNPcount = round((window/cov_sites)*nhet, 3)
	except ZeroDivisionError:
		SNPcount = 0.0
	output_f = Path(path, filename + '.heterozygosity.txt')
	with open(output_f, 'a') as out:
		out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome, start, end, cov_sites, nhet, SNPcount))




if __name__ == "__main__":
	args = parser.parse_args()
	filename = Path(args.vcf).stem.replace('.vcf', '')
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	sequences = sequences_bam(args.bam)
	max_depth = average_genome_coverage(args.cov)
	binned_heterozygosity = calculate_binned_heterozygosity(sequences, args.w, args.dp, max_depth, args.bam, args.vcf, args.gq, args.q, filename, args.o)
