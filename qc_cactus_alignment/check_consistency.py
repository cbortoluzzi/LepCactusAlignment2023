#!/usr/bin/env python


# Author : @cb46


import argparse
from Bio import AlignIO
from collections import defaultdict


parser = argparse.ArgumentParser(description = 'Evaluate consistency of single copy busco genes or single copy orthogroups')
parser.add_argument('--maf', help = 'Alignment in multiple alignment format')
parser.add_argument('--bed', help = 'A tab delimited file with genomic coordinate of single copy busco genes or single copy orthogroups')
parser.add_argument('--bp', help = 'Number of base pairs to take upstream and downstrean the genomic coordinate of a gene', type=int)
parser.add_argument('--refGenome', help = 'Reference genome')


def species_single_copy_coordinates(bed, base_pair, refGenome):
	mydict = {}
	with open(bed) as f:
		for line in f:
			contig, start, end, length, species, genome = line.strip().split()
			start, end = int(start), int(end)
			if base_pair:
				if genome == refGenome:
					mydict[genome] = [contig, start, end]
				else:
					new_start = start - base_pair
					new_end = end + base_pair
					mydict[genome] = [contig, new_start, new_end]
			else:
				mydict[genome] = [contig, start, end]
	return mydict



def parse_alignment(maf):
	mymaf = defaultdict(lambda: defaultdict(list))
	for multiple_alignment in AlignIO.parse(maf, "maf"):
		for seqrec in multiple_alignment:
			record, start, length, seq_size, strand, sequence = seqrec.id, seqrec.annotations['start'], seqrec.annotations['size'], seqrec.annotations['srcSize'], seqrec.annotations['strand'], str(seqrec.seq)
			species, contig = record.split('.', 1)
			contig = contig.lstrip('0')
			end = int(start) + int(length)
			seq_strand = ['+' if strand == 1  else '-']
			if seq_strand[0] == "+":
				new_start, new_end = int(start), int(end)
			else:
				new_start = int(seq_size) - end
				new_end = int(seq_size) - int(start)
			mymaf[species][contig].append((new_start, new_end, sequence))
	return mymaf



def check_consistency_alignment(mydict, mymaf, output_f):
	passed_consistency_check = defaultdict(lambda: defaultdict(list))
	to_be_checked = defaultdict(lambda: defaultdict(list))
	for species, (contig, start, end) in mydict.items():
		for interval in mymaf[species][contig]:
			begin, stop, seq = interval[0], interval[1], interval[2]
			# The alignment is found
			if begin in range(start, end + 1) and stop in range(start, end + 1):
				passed_consistency_check[species][contig, start, end].append('complete')
			# Only the end of the alignment is found
			elif stop in range(start, end + 1):
				to_be_checked[species][contig, start, end].append('partial')
			# Only the begin of the alignment is found
			elif begin in range(start, end + 1):
				to_be_checked[species][contig, start, end].append('partial')
	return passed_consistency_check, to_be_checked



def save_consistency_check_to_file(passed_consistency_check, to_be_checked, mydict):
	with open(output_f, 'w') as out_f:
		for species_name in passed_consistency_check:
			for region in passed_consistency_check[species_name]:
				type_region = passed_consistency_check[species_name][region][0]
				out_f.write('{}\t{}\t{}\t{}\t{}\n'.format(species_name, region[0], region[1], region[2], type_region))

		for key in to_be_checked:
			if key not in passed_consistency_check:
				for value in to_be_checked[key]:
					type_region = to_be_checked[key][value][0]
					out_f.write('{}\t{}\t{}\t{}\t{}\n'.format(key, value[0], value[1], value[2], type_region))

		for name in mydict.keys():
			if name not in passed_consistency_check and name not in to_be_checked:
				out_f.write('{}\n'.format('###'))
				out_f.write('{}\t{}\n'.format(name, 'Absent'))



if __name__ == "__main__":
	args = parser.parse_args()
	species_coordinates = species_single_copy_coordinates(args.bed, args.bp, args.refGenome)
	maf_f = parse_alignment(args.maf)
	output_f = args.maf.replace('.maf', '.maf.qc')
	passed_consistency_check, to_be_checked = check_consistency_alignment(species_coordinates, maf_f, output_f)
	save_to_file = save_consistency_check_to_file(passed_consistency_check, to_be_checked, species_coordinates)


