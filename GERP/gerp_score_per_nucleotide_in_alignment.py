#!/usr/bin/env python


# Author: @cb46



import argparse
from Bio import AlignIO
from pathlib import Path
from collections import defaultdict


parser = argparse.ArgumentParser(description = 'Assign the neutral rate and rejected substitution score (i.e. GERP score) to each position in the alignment')
parser.add_argument('--maf', help = 'Alignment in multiple alignment (MAF) format')



def get_alignment_positions(maf, refGenome):
	list_starts, list_ends = [], []
	list_seq = {}
	myrange = defaultdict(list)
	for multiple_alignment in AlignIO.parse(maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(refGenome):
				srcChrom, srcStart, srcSize, srcSeq = seqrec.id.split('.')[1], seqrec.annotations['start'], seqrec.annotations['size'], list(seqrec.seq)
				srcEnd = srcSize + srcStart
				count = -1
				for position in range(srcStart, srcEnd):
					count += 1
					nucleotide = srcSeq[count]
					myrange[srcChrom].append([position, nucleotide])
				list_starts.append(srcStart)
				list_ends.append(srcEnd)
	# Let's take into account the missing blocks because GERP++ outputs a value also for these blocks
	missing_blocks = list(zip(list_ends, list_starts[1:]))
	for (start, end) in missing_blocks:
		for i in range(start, end):
			# Since we do not have an alignment, we will replace the nucleotide with a gap (i.e. '-')
			myrange[srcChrom].append([i, '-'])
	return myrange




def gerp_score(gerp):
	mygerp = []
	with open(gerp) as gerp_f:
		for line in gerp_f:
			neutral_rate, rejected_substitution_score = map(float, line.strip().split())
			mygerp.append((neutral_rate, rejected_substitution_score))
	return mygerp




def assign_gerp_score_to_positions(mygerp, myrange, output_file):
	with open(output_file, 'w') as f:
		for chromosome in myrange:
			# Make sure that all positions are sorted !!
			list_positions = sorted(myrange[chromosome])
			gerp_score_position = list(zip(list_positions, mygerp))
			for element in gerp_score_position:
				begin, nucleotide = element[0][0], element[0][1]
				neutral_rate, rejected_substitution_score = element[1][0], element[1][1]
				stop = begin + 1
				# Let's remove positions for which we do not have a nucleotide (these are the gaps between blocks in the multiple sequence alignment)
				if nucleotide != '-':
					f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome, begin, stop, nucleotide, neutral_rate, rejected_substitution_score))




if __name__ == "__main__":
	args = parser.parse_args()
	refGenome = Path(args.maf).name.split('.')[0]
	# GERP++ output
	gerp = str(Path(args.maf))+".4d.rates"
	output_file = gerp + '.bed'
	list_positions = get_alignment_positions(args.maf, refGenome)
	neutral_rate_rejected_substitution_score = gerp_score(gerp)
	assigned_gerp_score = assign_gerp_score_to_positions(neutral_rate_rejected_substitution_score, list_positions, output_file)

  
