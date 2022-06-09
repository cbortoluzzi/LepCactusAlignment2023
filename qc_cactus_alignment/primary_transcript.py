#!/usr/bin/env python


# Author : @cb46


import argparse
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Select longest transcript from a multi fasta file')
parser.add_argument('--pep', help = 'Proteome in FASTA format')


def select_primary_transcript(pep_f):
	# Set output directory and output file
	path = Path(pep_f).parents[0]
	file_name = Path(pep_f).name
	# Create folders in case they don't exist
	# One folder for the longest transcripts
	longest_transcript_output_directory = Path(path, 'primary_transcripts')
	longest_transcript_output_directory.mkdir(parents=True, exist_ok=True)
	# One folder for the genomic coordinates of the longest transcript (these files will be used later on in the quality check step)
	genomic_coordinates_longest_transcript = Path(path, 'primary_transcripts_coordinates')
	genomic_coordinates_longest_transcript.mkdir(parents=True, exist_ok=True)
	longest_transcript = defaultdict(list)
	# Parse proteome in FASTA format
	for record in SeqIO.parse(pep_f, "fasta"):
		gene = record.description.split()[3].split(':')[1]
		coordinates = record.description.split()[2]
		tol_id, chromosome, start, end, strand = coordinates.split(':')
		len = int(end) - int(start)
		seq = str(record.seq)
		longest_transcript[gene].append([chromosome, start, end, len, seq])
	with open(Path(longest_transcript_output_directory, file_name), 'w') as fasta, open(Path(genomic_coordinates_longest_transcript, file_name.replace('.pep.fa', '.bed')), 'w') as bed:
		for gene_name in longest_transcript:
			# Select longest transcript
			longest = max(longest_transcript[gene_name], key=lambda a: a[3])
			fasta.write('{}\n{}\n'.format('>' + gene_name, longest[-1]))
			bed.write('{}\t{}\t{}\t{}\n'.format(gene_name, longest[0], longest[1], longest[2]))



if __name__ == "__main__":
	args = parser.parse_args()
	primary_transcript = select_primary_transcript(args.pep)

  
  
