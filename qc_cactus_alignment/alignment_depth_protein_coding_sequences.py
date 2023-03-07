#!/usr/bin/env python



# Author : @cb46



import argparse
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Assign alignment depth to each unique, non-overlapping protein-coding sequence')
parser.add_argument('--bed', help = 'Genomic coordinates of protein-coding sequences in BED format')
parser.add_argument('--wig', help = 'Path to alignment depth of protein-coding sequences in WIG format')
parser.add_argument('--g', help = 'Name of reference genome')




def unique_protein_coding_sequences(cds_f):
	mydict = {}
	with open(cds_f) as f:
		for line in f:
			chromosome, start, end = line.strip().split()
			start = int(start)
			end = int(end)
			length = end - start
			mydict[chromosome, start, end] = length
	return mydict



def alignment_depth_cds(wig_f, mydict, genome, path):
	wigD = defaultdict(list)
	cds_depth = defaultdict(int)
	for wig in wig_f:
		with open(wig) as f:
			n = 0
			for line in f:
				if line.startswith('fixedStep'):
					fixedStep, chromosome, start, step = line.strip().split()
					chrom = chromosome.split('=')[1]
					begin = start.split('=')[1]
					stepsize = step.split('=')[1]
					begin = int(begin)
					stepsize = int(stepsize)
				else:
					n += stepsize
					depth = int(line.strip())
					wigD[chrom].append([n, depth])

	for (chromosome, start, end) in mydict:
		for (position, depth) in wigD[chromosome]:
			if position >= start and position <= end:
        # Retain only positions with an alignment depth above 0
				if depth > 0:
					cds_depth[chromosome, start, end] += 1

	output = Path(path, genome + '.alignmentDepth.txt')
	with open(output, 'w') as outFile:
		for key in cds_depth:
			outFile.write('{}\t{}\t{}\t{}\n'.format(key[0], key[1], key[2], cds_depth[key]))



if __name__ == "__main__":
	args = parser.parse_args()
	wig_f = list(Path(args.wig).rglob('*.wig'))
	cds = unique_protein_coding_sequences(args.bed)
	path = Path(args.bed).parents[0]
	alignment_depth = alignment_depth_cds(wig_f, cds, args.g, path)


	
