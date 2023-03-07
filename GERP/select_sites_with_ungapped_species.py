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
parser.add_argument('--n', help = 'Maximum number of gaps per site [default = 3]', type = int, default = 3)
parser.add_argument('--nsp', help = 'Minimum number of species per block [default = 3]', type = int, default = 3)
parser.add_argument('--o', help = 'Output directory')



def multiple_sequence_alignment(input):
	mymaf = defaultdict(list)
	maf, species_name = input[0], input[1]
	for multiple_alignment in AlignIO.parse(maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(species_name):
				srcStart = seqrec.annotations["start"]
				srcSeq = str(seqrec.seq)
				mymaf[species_name, srcStart, srcSeq]
			else:
				tgtId = seqrec.id
				tgtSeq = str(seqrec.seq)
				mymaf[species_name, srcStart, srcSeq].append([tgtId, tgtSeq])
	return mymaf



def ungapped_sequences(mymaf, max_num_gap, chromosome, species_name, min_num_species, output_file):
	list_pos = []
	for key in mymaf[0]:
		targetSpecies = [i[0] for i in mymaf[0][key]]
		# We want to retain only alignment blocks with at least 3 species
		if len(targetSpecies) >= min_num_species:
			list_sequences = [i[-1] for i in mymaf[0][key]]
			sequences = list(zip(*list_sequences))
			n = 0
			start = key[1]
			for seq in sequences:
				gapCount = seq.count('-')
				nuclCount = len(seq) - gapCount
				n += 1
				position = start + n
				# Within each alignment block, we want to retain only sites with less than N gaps (N = 3 by default)
				if gapCount < max_num_gap:
					list_pos.append(position)

	with open(output_file, 'w') as out:
		groups = groupby(list_pos, key=lambda item, c=count():item-next(c))
		tmp = [list(g) for k, g in groups]
		for elem in tmp:
			begin = elem[0]
			stop = elem[-1]
			out.write('{}\t{}\t{}\n'.format(chromosome, begin, stop))




if __name__ == "__main__":
	args = parser.parse_args()
	# Generate folder if it doesn't exist
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	species_name, chromosome = Path(args.maf).stem.split('.')
	# Start multiprocessing ...
	num_workers = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(num_workers)
	input = list(zip([args.maf], [species_name]))
	parse_maf = pool.map(multiple_sequence_alignment, input)
	pool.close()
	pool.join()
	# Define output file
	output_name = Path(args.maf).stem + '.bed'
	output_file = Path(args.o, output_name)
	# Select sites with no more than 3 gaps
	sequence_range = ungapped_sequences(parse_maf, args.n, chromosome, species_name, args.nsp, output_file)

	
	
