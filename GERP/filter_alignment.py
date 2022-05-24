#!/usr/bin/env python


# Author : @cb46



import time
import argparse
import multiprocessing
from Bio import AlignIO
from pathlib import Path
from collections import defaultdict, Counter




parser = argparse.ArgumentParser(description = 'Prepare multiple sequence alignment for GERP++ analysis')
parser.add_argument('--maf', help = 'Alignment in multiple aligment format (MAF)')
parser.add_argument('--min_sp', help = 'Minimum number of species to have in each alignment block [default = 3]', type = int, default = 3)
parser.add_argument('--max_gap', help = 'Maximum number of species with gaps in sequence [default = 3]', type = int, default = 3)
parser.add_argument('--o', help = 'Output directory')




def multiple_sequence_alignment(input):
	maf, refGenome = input[0], input[1]
	mymaf = defaultdict(list)
	for multiple_alignment in AlignIO.parse(maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(refGenome):
				srcChrom, srcStart, srcStrand, srcSize, srcLength, srcSeq = seqrec.id, seqrec.annotations["start"], seqrec.annotations["strand"], seqrec.annotations["srcSize"], seqrec.annotations["size"], str(seqrec.seq)
				src_strand = ['+' if srcStrand == 1  else '-']
				mymaf[srcChrom, srcStart, srcLength, src_strand[0], srcSize, srcSeq]
			else:
				tgtChrom, tgtStart, tgtStrand, tgtSize, tgtLength, tgtSeq = seqrec.id, seqrec.annotations["start"], seqrec.annotations["strand"], seqrec.annotations["srcSize"], seqrec.annotations['size'], str(seqrec.seq)
				tgt_strand = ['+' if tgtStrand == 1  else '-']
				mymaf[srcChrom, srcStart, srcLength, src_strand[0], srcSize, srcSeq].append((tgtChrom, tgtStart, tgtLength, tgt_strand[0], tgtSize, tgtSeq))
	return mymaf



def filter_gap_sequences(mymaf, output_maf, min_sp, max_gap):
	with open(output_maf, 'w') as out:
		for query in mymaf[0]:
			if mymaf[0][query]:
				list_species = [i[0].split('.')[0] for i in mymaf[0][query]]
				number_unique_species = Counter(list_species)
				if len(number_unique_species) >= min_sp:
					list_seq = [i[-1] for i in mymaf[0][query]]
					count_gap = [i.count('-') for i in list_seq]
					seq_to_retain = [i for i in count_gap if i > 0]
					# Filter out alignment blocks with less than 3 species and with more than 3 gapped sequences
					if len(seq_to_retain) < max_gap:
						out.write('{}\n'.format('a'))
						query_maf = '\t'.join(map(str, query))
						out.write('{}\t{}\n'.format('s', query_maf))
						for target in mymaf[0][query]:
							target_maf = '\t'.join(map(str, target))
							out.write('{}\t{}\n'.format('s', target_maf))
					out.write('\n')




if __name__ == "__main__":
	args = parser.parse_args()
	# Check if directory exists, otherwise create it
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	# Obtain species name from alignment file
	sample = Path(args.maf).stem
	# Start multiprocessing ...
	num_workers = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(num_workers)
	input = list(zip([args.maf], [sample]))
	maf_to_dict = pool.map(multiple_sequence_alignment, input)
	pool.close()
	pool.join()
	# Filter multiple sequence alignment
	output_maf = Path(path, sample + '.maf.gerp')
	ungapped_maf = filter_gap_sequences(maf_to_dict, output_maf, args.min_sp, args.max_gap)

