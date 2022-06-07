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
parser.add_argument('--score', help = 'Minimum percentage identity of candidate conserved element', type = int, default = 70)
parser.add_argument('--length', help = 'Minimum length of candidate conserved element', type = int, default = 50)



def parse_multiple_sequence_alignment(maf, refGenome):
	mymaf = defaultdict(list)
	for multiple_alignment in AlignIO.parse(maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(refGenome):
				srcChrom = seqrec.id.split('.')[1]
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
						mymaf[srcChrom, srcStart, srcLen, srcStrand, srcSize, srcSeq].append([targetChromosome, tgtStart, tgtLen, tgtStrand, tgtSize, tgtSeq])
				except ValueError:
					pass
	return mymaf



def ungapped_sequences(mymaf):
	stepsize = 1
	ungapped_seq = defaultdict(list)
	for query in mymaf:
		for target in mymaf[query]:
			query_seq = query[-1]
			target_seq = target[-1]
			query_seq_array = np.array(list(query_seq))
			target_seq_array = np.array(list(target_seq))
			# Check that the two sequences have the same length
			assert len(query_seq_array) == len(target_seq_array), "Query and target sequence have different length"
			# Define index for query sequence
			index_q = [i for i, idx in enumerate(query_seq)]
			index_q_array = np.array(list(index_q))
			# Define index for target sequence: we do not add 1 to the index when in presence of gaps, on the contrary we add -1
			i = 0
			index_t = []
			for nucleotide in target_seq:
				if nucleotide != '-':
					index_t.append(i)
					i += 1
				else:
					index_t.append(-1)
			index_t_array = np.array(list(index_t))
			zipped_list = list(zip(index_q_array, index_t_array, query_seq_array, target_seq_array))
			# Get consecutive indeces in target sequence
			consecutive = np.split(index_t_array, np.where(np.diff(index_t_array) != stepsize)[0]+1)
			# Let's remove the -1 that indicates gaps
			for position in consecutive:
				new_position = [x for x in position if x != -1]
				if new_position:
					start, end = new_position[0], new_position[-1]
					# Get start and end position of query sequence, nucleotide sequence of query and nucleotide sequence of target
					list_index_q, list_nucl_q, list_nucl_t = [], [], []
					for (idx_q, idx_t, nucl_q, nucl_t) in zipped_list:
						if idx_t in range(start, end + 1) and idx_t != -1:
							list_index_q.append(idx_q)
							list_nucl_q.append(nucl_q)
							list_nucl_t.append(nucl_t)
					ungapped_query_seq = ''.join(list_nucl_q)
					ungapped_target_seq = ''.join(list_nucl_t)
					# Redefine start and end position for query ...
					start_q = query[1] + list_index_q[0]
					end_q = query[1] + list_index_q[-1]
					# Redefine start and end position for target ...
					start_t = target[1] + start
					end_t = target[1] + end
					ungapped_seq[query[0], start_q, end_q + 1, query[3], ungapped_query_seq].append([target[0], start_t, end_t + 1, target[3], target[4], ungapped_target_seq])
	return ungapped_seq



def filtering_ungapped_sequence(ungapped_seq, min_length, min_score):
	myelem = defaultdict(list)
	for query in ungapped_seq:
		for target in ungapped_seq[query]:
			query_seq = query[-1]
			target_seq = target[-1]
			len_query_seq = query[2] - query[1]
			len_target_seq = target[2] - target[1]
			# Check that two sequences have the same length
			assert len_query_seq == len_target_seq, "Sequences have different length"
			if len_query_seq >= min_length:
				# Change coordinates based on strand
				strand = ['+' if target[3] == 1  else '-']
				if strand[0] == "+":
					start = target[1]
					end = target[2]
				else:
					start = target[4] - target[2]
					end = target[4] - target[1]
				score = []
				percentage_identity = calculate_percentage_identity(query_seq, target_seq, score)
				# Retain pair of sequences with a minimum identity score of N
				if percentage_identity >= min_score:
					# Percentage of repeats in query and target sequence
					fraction_repeats_q = percentage_sequence_repeats(query_seq)
					fraction_repeats_t = percentage_sequence_repeats(target_seq)
					myelem[query[0], query[1], query[2], len_query_seq, '+', fraction_repeats_q, query_seq].append([target[0], start, end, strand[0], fraction_repeats_t, percentage_identity, target_seq])
	return myelem



def calculate_percentage_identity(query_seq, target_seq, score):
	len_t_seq = len(target_seq)
	for (nuclA, nuclB) in list(zip(query_seq, target_seq)):
		count = 0
		if nuclA.upper() == nuclB.upper():
			count += 1
			score.append(count)
	sum_score = sum(score)
	perc_identity = round((sum_score / len_t_seq) * 100, 2)
	return perc_identity



def percentage_sequence_repeats(sequence):
	len_seq = len(sequence)
	count_repeats = sum(1 for c in sequence if c.islower())
	fraction_repeats = round((count_repeats / len_seq) * 100, 2)
	return fraction_repeats



def write_ancestral_conserved_elements(myelem, refGenome, path):
	output_file = Path(path, 'ancestral_conserved_elements_' + refGenome + '.txt')
	with open(output_file, 'w') as output:
		for query in myelem:
			query_v = '\t'.join(map(str, query))
			for target in myelem[query]:
				target_v = '\t'.join(map(str, target))
				output.write('{}\t{}\n'.format(target_v, query_v))


				
				

if __name__ == "__main__":
	args = parser.parse_args()
	# Create directory based reference sequence name
	path = Path(args.refSequence)
	path.mkdir(parents = True, exist_ok = True)
	# Obtain relevant information from multiple sequence alignment
	maf_f = parse_multiple_sequence_alignment(args.maf, args.refGenome)
	# Obtain ungapped sequences for the target species
	ungapped_seq = ungapped_sequences(maf_f)
	# Filter ungapped sequences
	filtered_seq = filtering_ungapped_sequence(ungapped_seq, args.length, args.score)
	save_to_file = write_ancestral_conserved_elements(filtered_seq, args.refGenome, path)

