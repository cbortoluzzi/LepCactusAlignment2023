#!/usr/bin/env python



# Author : @cb46


import argparse
import multiprocessing
from Bio import AlignIO
from pathlib import Path
from collections import defaultdict
from itertools import groupby, count



parser = argparse.ArgumentParser(description = 'Identify sites in the alignment with at least 3 ungapped sequences')
parser.add_argument('--maf', help = 'Input maf file')
parser.add_argument('--o', help = 'Output directory')



def multiple_sequence_alignment(input):
	mymaf = defaultdict(list)
	maf, species_name = input[0], input[1]
	for multiple_alignment in AlignIO.parse(maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(species_name):
				srcStart, srcSeq = seqrec.annotations["start"], str(seqrec.seq)
				mymaf[species_name, srcStart, srcSeq]
			else:
				tgtId, tgtSeq = seqrec.id, str(seqrec.seq)
				mymaf[species_name, srcStart, srcSeq].append([tgtId, tgtSeq])
	return mymaf


def ungapped_sequences(mymaf, species_name, output_file):
	list_pos = []
	for key in mymaf[0]:
		if len(mymaf[0][key]) >= 3:
			sequences = [i[-1] for i in mymaf[0][key]]
			all_seqs = list(zip(*sequences))
			n = -1
			start = key[1]
			for seq in all_seqs:
				count_gap = seq.count('-')
				count_nucl = len(seq) - count_gap
				n += 1
				position = start + n
				if count_nucl >= 3:
					list_pos.append(position)


	with open(output_file, 'w') as out:
		groups = groupby(list_pos, key=lambda item, c=count():item-next(c))
		tmp = [list(g) for k, g in groups]
		for elem in tmp:
			begin = elem[0]
			stop = elem[-1]
			length = stop - begin
			sample, chromosome = species_name.split('.')
			if length > 0:
				out.write('{}\t{}\t{}\t{}\n'.format(chromosome, begin, stop, length))


if __name__ == "__main__":
	args = parser.parse_args()
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	species_name = Path(args.maf).stem
	# Start multiprocessing ...
	num_workers = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(num_workers)
	input = list(zip([args.maf], [species_name]))
	parse_maf = pool.map(multiple_sequence_alignment, input)
	pool.close()
	pool.join()
	output_file = Path(args.o, species_name + '.bed')
	sequence_range = ungapped_sequences(parse_maf, species_name, output_file)


