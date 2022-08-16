#!/usr/bin/env python



# Author : @cb46


import os
import vcf
import argparse
import subprocess
from pathlib import Path




parser = argparse.ArgumentParser(description = 'Polarize and filter sites based on phred-quality score, read depth, and genotype quality')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--cov', help = 'Genome-wide coverage')
parser.add_argument('--hal', help = 'Input hal file')
parser.add_argument('--refGenome', help = 'Name of reference genome')
parser.add_argument('--q', help = 'Minimum phred-quality score [default = 15]', type = int, default = 15)
parser.add_argument('--dp', help = 'Minimum read depth to retain a variant [default = 6]', type = int, default = 6)
parser.add_argument('--gq', help = 'Minimum genotype quality to retain a variant [default = 20]' , type = int, default = 20)
parser.add_argument('--o', help = 'Output directory')




def average_genome_coverage(coverage):
	with open(coverage) as f:
		for line in f:
			line = line.split('=')
			avg_genome_coverage = line[1].replace(' ','')
			max_depth = 2 * float(avg_genome_coverage)
	return max_depth



def filter_and_polarize_vcf_file(vcf_f, min_qual, min_depth, max_depth, min_gq, hal, refGenome, ancestor, output_file):
	vcf_reader = vcf.Reader(filename=vcf_f)
	vcf_writer = vcf.Writer(open(output_file, 'w'), vcf_reader)
	for record in vcf_reader:
		# Let's consider only bi-allelic SNPs that pass all filtering criteria
		if not record.FILTER and record.is_snp:
			# Filter based on phred-quality score
			if record.QUAL >= min_qual:
				for call in record.samples:
					read_depth = call['DP']
					genotype_quality = call['GQ']
					# Filter based on read depth and genotype quality
					if read_depth >= min_depth and read_depth <= max_depth and genotype_quality >= min_gq:
						chromosome = record.CHROM
						# We need to substract 1 from the position, because the VCF file is 1-based
						start = record.POS - 1
						try:
							# Obtain ancestral allele
							command = 'hal2maf --refSequence %s --start %s --length 1 --refGenome %s --onlyOrthologs --targetGenomes %s %s stdout' %(chromosome, start, refGenome, ancestor, hal)
							cmd = subprocess.check_output(command, stderr=subprocess.STDOUT, shell = True).decode()
							if len(cmd.split('\n')) == 5:
								output_reference = cmd.split('\n')[1].split()[-1]
								output_ancestor = cmd.split('\n')[2].split()[-1]
								assert record.REF == output_reference.upper(), "Reference allele is different from that in alignment"
								record.add_info('AA', output_ancestor)
								vcf_writer.write_record(record)
							else:
								output_reference = cmd.split('\n')[1].split()[-1]
								assert record.REF == output_reference.upper(), "Reference allele is different from that in alignment"
								record.add_info('AA', '.')
								vcf_writer.write_record(record)
						except subprocess.CalledProcessError as e:
							output = e.output.decode()
							continue


if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	output_file_name = Path(args.vcf).stem.replace('.vcf', '.filtered.vcf')
	output_file = Path(path, output_file_name)
	max_depth = average_genome_coverage(args.cov)
	# Obtain name of ancestor of target species
	command = 'halStats --parent %s %s' %(args.refGenome, args.hal)
	ancestor = subprocess.check_output(command, shell = True).decode().strip('\n')
	filtered_vcf = filter_and_polarize_vcf_file(args.vcf, args.q, args.dp, max_depth, args.gq, args.hal, args.refGenome, ancestor, output_file)
	
	

