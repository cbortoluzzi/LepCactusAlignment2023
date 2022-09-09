#!/usr/bin/env python



# Author : @cb46


import argparse
from pathlib import Path
from statistics import mean
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Define lepidoptera specific ancestral conserved elements')
parser.add_argument('--a', help = 'Ancestral conserved elements')
parser.add_argument('--b', help = 'Blastn output file')
parser.add_argument('--o', help = 'Output directory')



def blastn_output(blastn):
	mydict = defaultdict(list)
	with open(blastn) as f:
		for line in f:
			qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split()
			pident = float(pident)
			chromosome, ranges = qseqid.split(':')
			start, end = ranges.split('-')
			start = int(start)
			end = int(end)
			## Retain only significant matches
			if float(evalue) <= 0.05:
				mydict[chromosome, start, end].append(pident)
	return mydict



def classify_lepidoptera_conserved_elements(mydict, ancestral_f, path, output_file):
	output = Path(path, output_file)
	with open(output, 'w') as out:
		with open(ancestral_f) as f:
			for line in f:
				chromosome, start, length, strand, chromosome_len, perc_identity, perc_repeats, sequence = line.strip().split()
				start = int(start)
				end = start + int(length)
				try:
					## Type II are conserved in lepidoptera and have rudimentary homologous sequences in outgroup species
					## but the sequence conservation is significantly low to be detected as homologous
					if mydict[chromosome, start, end]:
						avg_pident = mean(mydict[chromosome, start, end])
						if avg_pident < 70:
							out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome, start, length, strand, chromosome_len, perc_identity, perc_repeats, sequence, 'type_II'))
					## Type I are conserved in lepidoptera but have no homologous sequence in other vertebrate outgroups
					if not mydict[chromosome, start, end]:
						out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome, start, length, strand, chromosome_len, perc_identity, perc_repeats, sequence, 'type_I'))
				except KeyError:
					continue




if __name__ == "__main__":
	args = parser.parse_args()
	# Create output directory if it doesn't exist
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	output_file = Path(args.b).stem + '.typeI_II.bed'
	blastn = blastn_output(args.b)
	type_I_II_elements = classify_lepidoptera_conserved_elements(blastn, args.a, path, output_file)

