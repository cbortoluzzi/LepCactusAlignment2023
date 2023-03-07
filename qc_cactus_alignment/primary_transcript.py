#!/usr/bin/env python


# Author : @cb46


import gzip
import argparse
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Select primary transcript before running OrthoFinder')
parser.add_argument('--pep', help = 'Protein file in FASTA format')
parser.add_argument('--o', help = 'Output directory')


def primary_transcript(protein, path, output_pep, output_longest_transcript):
	mydict = defaultdict(list)
	output_file = Path(path, output_pep)
	output_transcript = Path(path, output_longest_transcript)
	with gzip.open(protein, 'rt') as prot:
		for record in SeqIO.parse(prot, "fasta"):
			description = record.description.split()
			coordinates = description[2].split(':')
			gene = description[3].split(':')[1]
			transcript = description[4].split(':')[1]
			sequence = str(record.seq)
			start = int(coordinates[-3])
			end = int(coordinates[-2])
			length = end - start
			mydict[gene].append([coordinates[-4], start, end, length, coordinates[-1], sequence, transcript])
	with open(output_file, "w") as output_1, open(output_transcript, 'w') as output_2:
		for key in mydict:
			list_len = [i[3] for i in mydict[key]]
			longest = max(list_len)
			index = list_len.index(longest)
			longest_transcript = mydict[key][index]
			coordinates_long_transcript = ':'.join(map(str, longest_transcript[0:5]))
			header_fasta = '>' + key + ':' + coordinates_long_transcript
			output_1.write('{}\n{}\n'.format(header_fasta, longest_transcript[5]))
			output_2.write('{}\t{}\n'.format(key, longest_transcript[-1]))




if __name__ == "__main__":
	args = parser.parse_args()
	# Generate directory if it doesn't exist
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	output_pep = Path(args.pep).stem
	output_longest_transcript = Path(args.pep).stem.replace('.fa',  '.txt')
	longest_transcript = primary_transcript(args.pep, path, output_pep, output_longest_transcript)


