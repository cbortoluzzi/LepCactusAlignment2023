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
parser.add_argument('--min_identity', help = 'Minimum percentage identity of candidate conserved element', type = int, default = 60)
parser.add_argument('--min_length', help = 'Minimum length of candidate conserved element', type = int, default = 50)



def parse_multiple_sequence_alignment(maf, refGenome):
	"""
	Read a multiple sequence alignment in MAF format
	Input:
		maf : multiple sequence alignment
		refGenome : name of reference genome
	Output:
		mymaf : dictionary with relevant information on query and target sequence
	"""
	mymaf = defaultdict(list)
	for multiple_alignment in AlignIO.parse(maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(refGenome):
				srcChrom, srcStart, srcStrand, srcSize, srcSeq = seqrec.id.split('.', 1)[1], seqrec.annotations["start"], seqrec.annotations["strand"], seqrec.annotations["srcSize"], str(seqrec.seq)
				mymaf[refGenome, srcChrom, srcStart, srcStrand, srcSize, srcSeq]
			else:
				tgtChrom, tgtStart, tgtStrand, tgtSize, tgtSeq = seqrec.id.split('.', 1)[1], seqrec.annotations["start"], seqrec.annotations["strand"], seqrec.annotations["srcSize"], str(seqrec.seq)
				targetGenome = seqrec.id.split('.')[0]
				mymaf[refGenome, srcChrom, srcStart, srcStrand, srcSize, srcSeq].append((targetGenome, tgtChrom, tgtStart, tgtStrand, tgtSize, tgtSeq))
	return mymaf



def ungapped_sequences(mymaf):
	"""
	Split target sequence at gaps in order to obtain ungapped sequences
	Input:
		mymaf : dictionary with relevant information on query and target sequence
	Output:
		mydict : dictionary of ungapped sequences
	"""
	stepsize = 1
	mydict = defaultdict(list)
	for query in mymaf:
		for target in mymaf[query]:
			query_seq, target_seq = query[-1], target[-1]
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
					start_q = query[2] + list_index_q[0]
					end_q = query[2] + list_index_q[-1]
					# Redefine start and end position for target ...
					start_t = target[2] + start
					end_t = target[2] + end
					mydict[query[1], start_q, end_q + 1, query[3], query[4], ungapped_query_seq].append((target[0], target[1], start_t, end_t + 1, target[3], target[4], ungapped_target_seq))
	return mydict



def select_single_and_multiple_copy_sequences(mydict):
	"""
	Define sequences as single-copy or multiple-copy based on the number of sequences the query aligns to
	We perform this step before applying any additional filtering on the length and identity score
	Input:
		mydict : dictionary of ungapped sequences
	Output:
		single_copy : dictionary with all ungapped single-copy sequneces
		multiple_copy : dictionary with all ungapped multiple-copy sequences
	"""
	single_copy, multiple_copy = {}, {}
	for query in mydict:
		if mydict[query]:
			# If the query aligns to only one target sequence ...
			if len(mydict[query]) == 1:
				single_copy[query] = mydict[query]
			# otherwise ...
			else:
				multiple_copy[query] = mydict[query]
	return single_copy, multiple_copy



def select_conserved_elements(single_copy, multiple_copy, min_identity, min_len, refSequence, targetGenome):
	"""
	Perform additional filtering on ancestral (single/multiple-copy) sequences and retain only those that meet our criteria of length and % identity
	Input:
		single_copy, multiple_copy : dictionary of candidate ancestral sequences
		min_identity : minimum value for identity score between query and target (default: 60%)
		min_len : minimum value for sequence length (default: 50 bp)
		refSequence : name of reference sequence
		targetGenome : name of target genome
	Output:
		ancestral_conserved_elements tab delimited file with relevant information on all ancestral conserved elements
	"""
	single_copy_conserved = filtering_step(single_copy, min_identity, min_len)
	multiple_copy_conserved = filtering_step(multiple_copy, min_identity, min_len)

	# Write candidate conserved elements to file
	output = Path(refSequence, 'ancestral_conserved_elements_' + targetGenome + '.txt')
	single_copy_to_file = write_conserved_elements(single_copy_conserved, output, 'single_copy')
	multiple_copy_to_file = write_conserved_elements(multiple_copy_conserved, output, 'multiple_copy')

	return single_copy_conserved, multiple_copy_conserved



def filtering_step(dictionary, minimum_identity, minimum_length):
	"""
	Filter ungapped (single/multiple-copy) sequences based on minimum length and minimum identity score
	Input:
		dictionary : dictionary with ungapped sequences (i.e. single-copy, multiple-copy)
		minimum_identity : minimum identity score (integer)
		minimum_length : minimum length of sequence (integer)
	Output:
		mycons : dictionary with filtered, ungapped conserved sequences
	"""
	mycons = defaultdict(list)
	for query in dictionary:
		for target in dictionary[query]:
			seq_q, seq_t = query[-1], target[-1]
			# Change coordinates based on strand
			strand = ['+' if target[4] == 1  else '-']
			if strand[0] == "+":
				start, end = target[2], target[3]
			else:
				start = target[5] - target[3]
				end = target[5] - target[2]
			clipped_seq_q, clipped_seq_t = seq_q, seq_t
			length_q = query[2] - query[1]
			length_t = end - start
			assert length_q == length_t == len(clipped_seq_q) == len(clipped_seq_t), "Sequences have different length"
			if length_t >= minimum_length:
				score = []
				for (nuclA, nuclB) in list(zip(clipped_seq_q, clipped_seq_t)):
					count = 0
					if nuclA.upper() == nuclB.upper():
						count += 1
						score.append(count)
				sum_score = sum(score)
				perc_identity = round((sum_score / length_t) * 100, 2)
				if perc_identity >= minimum_identity:
					# Fraction repeats in target ...
					count_repeats_t = sum(1 for c in clipped_seq_t if c.islower())
					frac_repeats_t = round((count_repeats_t / length_t) * 100, 2)
					# ... and query
					count_repeats_q = sum(1 for c in clipped_seq_q if c.islower())
					frac_repeats_q = round((count_repeats_q / length_q) * 100, 2)
					mycons[query[0], query[1], query[2], length_q, "+", frac_repeats_q, clipped_seq_q].append((target[1], start, end, strand[0], frac_repeats_t, perc_identity, clipped_seq_t))
	return mycons



def write_conserved_elements(dictionary, output, type_element):
	"""
	Write ancestral conserved elements to file

	Input:
		dictionary : list of filtered, ungapped ancestral conserved elements
		output : name of output file
	"""
	with open(output, 'a') as bed:
		for key in dictionary:
			keys = '\t'.join(map(str, key))
			bed.write('{}\n'.format('###'))
			for value in dictionary[key]:
				values = '\t'.join(map(str, value))
				bed.write('{}\t{}\t{}\n'.format(values, keys, type_element))




if __name__ == "__main__":
	args = parser.parse_args()
	# Create directory based reference sequence name
	path = Path(args.refSequence)
	path.mkdir(parents = True, exist_ok = True)
	maf_name = Path(args.maf).name
	species = maf_name.replace('.maf', '')
	maf_to_dict = parse_multiple_sequence_alignment(args.maf, args.refGenome)
	ungapped_seq = ungapped_sequences(maf_to_dict)
	single_copy, multiple_copy = select_single_and_multiple_copy_sequences(ungapped_seq)
	single_copy_conserved, multiple_copy_conserved = select_conserved_elements(single_copy, multiple_copy, args.min_identity, args.min_length, args.refSequence, species)
