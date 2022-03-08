#!/usr/bin/env python


# Author : @cb46



import argparse
from Bio import AlignIO
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Filter alignment for gaps')
parser.add_argument('--maf', help = 'Alignment in multiple aligment format (MAF)')
parser.add_argument('--refGenome', help = 'Reference genome')



def split_by_chromosome(maf, refGenome):
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



def filter_gap_sequences(mymaf, refGenome, output):
	with open(output, 'w') as out:
		for query in mymaf:
			if mymaf[query]:
				list_seq = [i[-1] for i in mymaf[query]]
				count_gap = [i.count('-') for i in list_seq]
				seq_to_retain = [i for i in count_gap if i > 0]
				# Filter out alignment blocks with less than 3 species and with more than 3 gapped sequences
				if len(list_seq) >= 3 and len(seq_to_retain) < 3:
					out.write('{}\n'.format('a'))
					query_maf = '\t'.join(map(str, query))
					out.write('{}\t{}\n'.format('s', query_maf))
					for target in mymaf[query]:
						target_maf = '\t'.join(map(str, target))
						out.write('{}\t{}\n'.format('s', target_maf))
					out.write('\n')


					

if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.maf).parents[0]
	filename = Path(args.maf).stem +'.maf.input.gerp'
	output_file = Path(path, filename)
	split_maf = split_by_chromosome(args.maf, args.refGenome)
	ungapped_maf = filter_gap_sequences(split_maf, args.refGenome, output_file)

  
