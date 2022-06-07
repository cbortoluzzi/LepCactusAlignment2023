#!/usr/bin/env python



# Author : @cb46



# Import modules
import argparse
import numpy as np
from Bio import AlignIO
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Identify ancestral conserved elements in a multiple sequence alignment')
parser.add_argument('--maf', help = 'Multiple alignment (MAF)')
parser.add_argument('--refGenome', help = 'Name of reference genome')
parser.add_argument('--refSequence', help = 'Name of reference sequence within reference genome')
parser.add_argument('--score', help = 'Minimum identity score of candidate conserved element [default = 70]', type = int, default = 70)
parser.add_argument('--length', help = 'Minimum length of candidate conserved element [default = 100]', type = int, default = 100)



def parse_multiple_sequence_alignment(maf, refGenome):
	mymaf = defaultdict(list)
	for multiple_alignment in AlignIO.parse(maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(refGenome):
				srcChrom = seqrec.id
				srcStart = seqrec.annotations["start"]
				srcLen = seqrec.annotations["size"]
				srcStrand = '+'
				srcSize = seqrec.annotations["srcSize"]
				srcSeq = str(seqrec.seq)
				mymaf[srcChrom, srcStart, srcLen, srcStrand, srcSize, srcSeq]
			else:
				targetGenome, targetChromosome = seqrec.id.split('.', 1)
				tgtStart = seqrec.annotations['start']
				tgtLen = seqrec.annotations["size"]
				tgtStrand = seqrec.annotations["strand"]
				tgtSize = seqrec.annotations["srcSize"]
				tgtSeq = str(seqrec.seq)
				# Let's retain only sex chromosomes (X, Y, W, or Z) and autosomes
				try:
					if targetChromosome == 'W' or targetChromosome == 'Z' or targetChromosome == 'X' or targetChromosome == 'Y' or isinstance(int(targetChromosome), int):
						tgtChrom = targetGenome + '.' + targetChromosome
						mymaf[srcChrom, srcStart, srcLen, srcStrand, srcSize, srcSeq].append([tgtChrom, tgtStart, tgtLen, tgtStrand, tgtSize, tgtSeq])
				except ValueError:
					pass
	return mymaf




def filter_sequences(mymaf, minimum_length, minimum_identity):
	mydict = defaultdict(list)
	for query in mymaf:
		# Retain only alignment blocks with a minimum length of N
		if query[2] >= minimum_length:
			for target in mymaf[query]:
				query_seq = query[-1]
				target_seq = target[-1]
				len_q_seq = len(query_seq)
				len_t_seq = len(target_seq)
				assert len_q_seq == len_t_seq, "Sequences have different length"
				score = []
				percentage_identity = calculate_percentage_identity(query_seq, target_seq, score)
				# Retain pair of sequences with a minimum identity score of N
				if percentage_identity >= minimum_identity:
					# Percentage of repeats in query and target sequence
					fraction_repeats_q = percentage_sequence_repeats(query_seq)
					fraction_repeats_t = percentage_sequence_repeats(target_seq)
					# Percentage of gaps in target sequence
					fraction_gaps_t = percentage_sequence_gaps(target_seq)
					query_info = '\t'.join(map(str, query))
					# Retain only sequences in the target species that have less than 10% gaps
					if fraction_gaps_t <= 10:
						# Change coordinates based on strand
						strand = ['+' if target[3] == 1  else '-']
						if strand[0] == '+':
							tgtStart = target[1]
						else:
							tgtEnd = target[1] + target[2]
							tgtStart = target[4] - tgtEnd
						mydict[query_info, fraction_repeats_q].append([target[0], tgtStart, target[2], strand[0], target[4], target_seq, fraction_repeats_t, percentage_identity])
	return mydict


def calculate_percentage_identity(query_seq, target_seq, score):
	len_t_seq = len(target_seq)
	for (nuclA, nuclB) in list(zip(query_seq, target_seq)):
		count = 0
		if nuclA.upper() == nuclB.upper():
			count += 1
			score.append(count)
	sum_score = sum(score)
	percentage_identity = round((sum_score / len_t_seq) * 100, 2)
	return percentage_identity


def percentage_sequence_repeats(sequence):
	len_seq = len(sequence)
	count_repeats = sum(1 for c in sequence if c.islower())
	fraction_repeats = round((count_repeats / len_seq) * 100, 2)
	return fraction_repeats


def percentage_sequence_gaps(sequence):
	len_seq = len(sequence)
	count_gaps = sum(1 for c in sequence if c == '-')
	fraction_gaps = round((count_gaps / len_seq) * 100, 2)
	return fraction_gaps


def generate_alignment(mydict, output_f, path):
	# Save to file only if dictionary is not empty
	if mydict:
		with open(output_f, 'w') as out:
			for query in mydict:
				out.write('{}\n'.format('a'))
				out.write('{}\t{}\n'.format('s', query[0]))
				for target in mydict[query]:
					target_info = '\t'.join(map(str, target[:6]))
					out.write('{}\t{}\n'.format('s', target_info))
					targetSpecies, targetChromosome = target[0].split('.')
					targetEnd = target[1] + target[2]
					genomic_coordinates_f = Path(path, 'ancestral_conserved_elements_' + targetSpecies + '.txt')
					with open(genomic_coordinates_f, 'a') as output:
						output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(targetChromosome, target[1], targetEnd, target[2], target[3], target[6], target[7]))
			out.write('\n')




if __name__ == "__main__":
	args = parser.parse_args()
	# Create directory based reference sequence name
	path = Path(args.refSequence)
	path.mkdir(parents = True, exist_ok = True)
	# Output file
	output_f = Path(path, 'ancestral_conserved_elements_' + args.refSequence + '.maf')
	# Obtain relevant information from multiple sequence alignment
	maf_to_dict = parse_multiple_sequence_alignment(args.maf, args.refGenome)
	# Filter out sequences based on minimim length and minimum identity score
	filtered_seq = filter_sequences(maf_to_dict, args.length, args.score)
	maf_alignment = generate_alignment(filtered_seq, output_f, path)

