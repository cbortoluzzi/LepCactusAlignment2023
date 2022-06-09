#!/usr/bin/env python


# Author : @cb46


import argparse
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Select longest transcript of a gene as primary transcript')
parser.add_argument('--pep', help = 'Proteome in FASTA format')


def select_primary_transcript(pep_f, file_name, path_1, path_2):
	longest_transcript = defaultdict(list)
	# Parse proteome in FASTA format
	for record in SeqIO.parse(pep_f, "fasta"):
		gene = record.description.split()[3].split(':')[1]
		coordinates = record.description.split()[2]
		tol_id, chromosome, start, end, strand = coordinates.split(':')
		len = int(end) - int(start)
		seq = str(record.seq)
		longest_transcript[gene].append([chromosome, start, end, len, seq])

	genomic_coordinates_f = file_name.replace('.pep.fa', '.bed')
	with open(Path(path_1, file_name), 'w') as fasta, open(Path(path_2, genomic_coordinates_f), 'w') as bed:
		for gene_name in longest_transcript:
			# Select longest transcript
			longest = max(longest_transcript[gene_name], key=lambda a: a[3])
			fasta.write('{}\n{}\n'.format('>' + gene_name, longest[-1]))
			bed.write('{}\t{}\t{}\t{}\n'.format(gene_name, longest[0], longest[1], longest[2]))



if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.pep).parents[0]
	file_name = Path(args.pep).name
	# Create folder where to print primary transcripts
	longest_transcript_output_directory = Path(path, 'primary_transcripts')
	longest_transcript_output_directory.mkdir(parents=True, exist_ok=True)
	# Create folder where to print the genomic coordinates of the longest transcript
	genomic_coordinates_longest_transcript = Path(path, 'primary_transcripts_coordinates')
	genomic_coordinates_longest_transcript.mkdir(parents=True, exist_ok=True)
	# Identify longest transcript
	primary_transcript = select_primary_transcript(args.pep, file_name, longest_transcript_output_directory, genomic_coordinates_longest_transcript)

  
