#!/usr/bin/env python



# Author : @cb46


import vcf
import argparse
from pathlib import Path



parser = argparse.ArgumentParser(description = 'Assign GERP score to polarized high-confident heterozygous sites')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--bed', help = 'Species folder with tab delimited bed files with GERP score information')
parser.add_argument('--o', help = 'Output directory')




def gerp_score_per_variant(bed_f):
	mydict = {}
	for bed in bed_f:
		with open(bed) as f:
			for line in f:
				chromosome, start, end, nucleotide, neutral_rate, rejected_substitution_score = line.strip().split()
				start, end = int(start), int(end)
				# Consider only sites with a GERP score > 5 as these represent a compromise between excluding as many as possible of tentatively neutral sites,
				# while including as many as possible of potentially deleterious sites
				if float(rejected_substitution_score) >= 5:
					mydict[chromosome, start, end] = [nucleotide, neutral_rate, rejected_substitution_score]
	return mydict




def assign_gerp_score_to_polarized_variants(mydict, vcf_f, output_file):
	vcf_reader = vcf.Reader(filename=vcf_f)
	with open(output_file, 'w') as out:
		for record in vcf_reader:
			for call in record.samples:
				# Retain only sites with information on the ancestral allele
				if record.INFO['AA'][0]:
					ancestral_allele = (record.INFO['AA'][0]).upper()
					alternative_allele = record.ALT[0]
					# Retain only those alleles for which the ancestral allele is either equal to reference or alternative allele
					if ancestral_allele == alternative_allele or ancestral_allele == record.REF:
						# We need to do this because the VCF file is 1-based
						start = record.POS - 1
						end = start + 1
						if (record.CHROM, start, end) in mydict.keys():
							key = mydict[record.CHROM, start, end]
							if record.REF == key[0] and call['GT'] == '0/1':
								out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(record.CHROM, start, end, record.REF, record.ALT[0], ancestral_allele, key[1], key[2]))









if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	bed_files = list(Path(args.bed).rglob('*.conserved.bed'))
	output_f = Path(args.vcf).stem + '.GERPscore.bed'
	output_file = Path(path, output_f)
	scored_variants = gerp_score_per_variant(bed_files)
	scored_polarized_variants = assign_gerp_score_to_polarized_variants(scored_variants, args.vcf, output_file)
  
  
