#!/usr/bin/env python



# Author : @cb46


import argparse
import subprocess
import multiprocessing
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(description = 'Obtain alignment depth for a genome')
parser.add_argument('--hal', help = 'Input hal file')
parser.add_argument('--refGenome', help = 'Name of reference genome to scan')
parser.add_argument('--refSequence', help = 'Name of reference sequence')
parser.add_argument('--start', help = 'Coordinate within reference genome (or sequence if specified) to start at [default = 0]', type = int, default = 0)
parser.add_argument('--length', help = 'Length of the reference genome (or sequence if specified). If set to 0, the entire thing is converted [default = 0]', type = int, default = 0)
parser.add_argument('--step', help = 'Step size [default = 1]', type = int, default = 1)
parser.add_argument('--depth', help = 'Minimum number of other unique genomes each base pair aligns to', type = int)
parser.add_argument('--plot', help = 'Option to plot the number of unique genomes each base aligns to [default = Yes]', type = str, default = 'Yes')
parser.add_argument('--o', help = 'Output directory')



def alignment_depth(hal, refGenome, sequence, start, length, step, min_depth, path):
	# Count the number of other unique genomes each base aligns to, excluding ancestral genomes
	mydepth = defaultdict(list)
	# If length == 0, then the entire genome is considered
	if length == 0:
		command = 'halAlignmentDepth --noAncestors --refSequence %s --step %d %s %s' %(sequence, step, hal, refGenome)
	# Otherwise, a region is considered
	else:
		command = 'halAlignmentDepth --noAncestors --refSequence %s --step %d --start %d --length %d %s %s' %(sequence, step, start, length, hal, refGenome)
	cmd = subprocess.check_output(command, shell = True).decode()
	outcmd = cmd.split('\n')
	for line in outcmd:
		if line:
			if line.startswith('fixedStep'):
				nbases, ncov = 0, 0
				fixedStep, chrom, start_pos, step_pos = line.strip().split()
				position = start_pos.split('=')[1]
				step_size = step_pos.split('=')[1]
				position, step_size = int(position), int(step_size)
			else:
				position += step_size
				depth = int(line)
				nbases += 1
				if depth >= min_depth:
					ncov += 1
					mydepth[sequence].append((position - 1, depth))
	# Print number of aligned bases
	output_f = Path(path, refGenome + '.tsv')
	with open(output_f, 'a') as tsv:
		end = start + nbases
		try:
			ncov = ncov
		except IndexError:
			ncov = 0
		tsv.write('{}\t{}\t{}\t{}\t{}\n'.format(sequence, start, end, nbases, ncov))
	return mydepth


def plot_halAlignmentDepth(mydepth, path, refGenome, sequence):
	x, y = [], []
	for key in mydepth[0]:
		x.append([i[0] for i in mydepth[0][key]])
		y.append([i[1] for i in mydepth[0][key]])
		# Plot depth versus position in the sequence
		plt.figure(figsize=(10,5))
		plt.plot(x[0], y[0])
		plt.xlabel('Position along the genome (bp)')
		plt.ylabel('Number of other unique genomes')
		plt.ylim(0, 89)
		fig = Path(path, refGenome + '.' + sequence + '.pdf')
		plt.savefig(fig, dpi = 500, bbox_inches = 'tight')




if __name__ == "__main__":
	args = parser.parse_args()
	# Create output directory if it doesn't exist
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	# We will run the halAlignmentDepth utility using multiprocessing
	num_workers = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(num_workers)
	wiggle = pool.starmap(alignment_depth, [(args.hal, args.refGenome, args.refSequence, args.start, args.length, args.step, args.depth, path)])
	pool.close()
	pool.join()
	# If --plot is chosen, then plot alignment depth
	if args.plot == 'Yes':
		plot_alignmentdepth = plot_halAlignmentDepth(wiggle, path, args.refGenome, args.refSequence)

