#!/usr/bin/env python



# Author : @cb46


import vcf
import argparse
import subprocess
from pathlib import Path



parser = argparse.ArgumentParser(description = 'Filter variants based on phred-quality score, read depth, and genotype quality')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--cov', help = 'Genome-wide coverage as estimated by samtools depth')
parser.add_argument('--q', help = 'Minimum phred-quality score [default = 15]', type = int, default = 15)
parser.add_argument('--dp', help = 'Minimum read depth to retain a variant [default = 6]', type = int, default = 6)
parser.add_argument('--gq', help = 'Minimum genotype quality to retain a variant [default = 20]' , type = int, default = 20)
parser.add_argument('--o', help = 'Output directory')




def average_genome_coverage(coverage):
	with open(coverage) as f:
		for line in f:
			line = line.split('=')
			avg_genome_coverage = line[1].replace(' ','')
			# Our maximum coverage is 2 times the average genome-wide coverage
			max_depth = 2 * float(avg_genome_coverage)
	return max_depth



def filter_vcf_file(vcf_f, min_qual, min_depth, max_depth, min_gq, output_file):
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
						vcf_writer.write_record(record)



if __name__ == "__main__":
	args = parser.parse_args()
	# Generate directory if it doesn't exist
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	# Set output file
	output_file_name = Path(args.vcf).stem.replace('.vcf', '.filtered.vcf')
	output_file = Path(path, output_file_name)
	max_depth = average_genome_coverage(args.cov)
	# Filter VCF file
	filter_vcf = filter_vcf_file(args.vcf, args.q, args.dp, max_depth, args.gq, output_file)
	
