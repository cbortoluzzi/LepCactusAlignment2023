#!/usr/bin/env python


# Author : @cb46


import gzip
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description = 'Change the header of a FASTA file')
parser.add_argument('--fa', help = 'FASTA file')
parser.add_argument('--report', help = 'Assembly report')
parser.add_argument('--out', help = 'Output FASTA file')



def parse_assembly_report(report):
	mycontig = {}
	with open(report) as f:
		for line in f:
			if not line.startswith('#'):
				line = line.strip().split()
				if line[1] == "assembled-molecule":
					chromosome_num, ncbi_contig = line[2], line[4]
					mycontig[ncbi_contig] = chromosome_num
				else:
					ncbi_scaffold = line[4]
					mycontig[ncbi_scaffold] = ncbi_scaffold
	return mycontig



def change_header_fasta(mycontig, fasta, output):
	with open(fasta) as handle, open(output, 'w') as out:
		records = SeqIO.parse(handle, 'fasta')
		for record in records:
			record.id = mycontig[record.id]
			record.description = ''
			SeqIO.write(record, out, 'fasta')


if __name__ == "__main__":
	args = parser.parse_args()
	assembly_report = parse_assembly_report(args.report)
	change_header = change_header_fasta(assembly_report, args.fa, args.out)

  
