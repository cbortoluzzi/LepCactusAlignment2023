#!/usr/bin/env python


# Author : @cb46


import os
import argparse
from Bio import SeqIO
import subprocess
from ete3 import Tree
from pathlib import Path
from Bio import AlignIO
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Check consistency of Orthogroups genomic coordinates with cactus alignment')
parser.add_argument('--refGenome', help = 'Name of species to use as query')
parser.add_argument('--t', help = 'Phylogenetic tree used as guide tree in cactus')
parser.add_argument('--f', help = 'Tab delimited species list file')
parser.add_argument('--ort', help = 'List of single copy orthogroups, one per line')
parser.add_argument('--hal', help = 'Input hal file')
parser.add_argument('--o', help = 'Output directory')



main_path = '/lustre/scratch123/tol/teams/durbin/users/cb46/whole_genome_alignment/lepidoptera/89way/202201/quality_check/orthofinder/proteomes/primary_transcripts'


def pairwise_comparisons(refGenome, tree, path):
	list_nodes = []
	tree = Tree(tree, format = 1)
	# Re-root the guide tree to the outgroup
	tree.set_outgroup('tinea_trinotella_gca905220615v1')
	for node in tree.traverse('postorder'):
		if node.is_leaf():
			if node.name != refGenome:
				output_directory = Path(path, refGenome + '_vs_' + node.name)
				output_directory.mkdir(parents=True, exist_ok=True)
				list_nodes.append(node.name)
	mypair = list(map(lambda e: (e, refGenome), list_nodes))
	return mypair



def get_species_name_tol_id(species_list):
	species_d = {}
	with open(species_list) as f:
		for line in f:
			tol_id, pclass, order, superfamily, family, latin_name, assembly = line.strip().split()
			genome = latin_name.lower() + '_' + assembly.replace('GCA_', 'gca').replace('.', 'v')
			species_d[genome] = [tol_id, assembly]
	return species_d




def get_single_copy_orthogroups(orthogroups, query, target, main_path, species_d, hal, path):
	with open(orthogroups) as f:
		for line in f:
			orthogroup_id = line.strip()
			gene_names = get_gene_name(orthogroup_id, query, target, main_path)
			for key, value in gene_names.items():
				query = key[1]
				target = key[2]
				output_maf = Path(path, query + '_vs_' + target, key[0] + '.maf')
				output_qc = Path(path, query + '_vs_' + target, key[0] + '.qc')
				# Query
				qSeq, qStart, qEnd, qLength = value[0], value[1], value[2], value[3]
				qGenome = query
				# Target species
				tSeq, tStart, tEnd = value[4], value[5], value[6]
				tGenome = target
				maf = pairwise_alignment(qSeq, qGenome, qStart, qEnd, qLength, tGenome, hal, output_maf, output_qc, tSeq, tStart, tEnd)
	return gene_names



def get_gene_name(orthogroup_id, query, target, main_path):
	dictionary = {}
	pep_query_vs_pep_target = Path(main_path, 'OrthoFinder/Results_Oct31/Orthologues/', 'Orthologues_' + query + '.pep', query + '.pep__v__' + target + '.pep.tsv')
	if Path(pep_query_vs_pep_target).is_file():
		with open(pep_query_vs_pep_target) as f:
			next(f)
			for line in f:
				line = line.strip().split()
				if orthogroup_id == line[0]:
					gene_query, qSeq, qStart, qEnd, qLength, qStrand = line[1].split('_')
					gene_target, tSeq, tStart, tEnd, tLength, tStrand = line[2].split('_')
					dictionary[orthogroup_id, query, target] = [qSeq, qStart, qEnd, qLength, tSeq, tStart, tEnd]
	return dictionary



def pairwise_alignment(qSeq, qGenome, qStart, qEnd, qLength, tGenome, hal, output_maf, output_qc, tSeq, tStart, tEnd):
	mypos = defaultdict(list)
	cmd = 'hal2maf --refSequence %s --refGenome %s --start %s --length %s --noAncestors --onlyOrthologs --targetGenomes %s %s %s' %(qSeq, qGenome, qStart, qLength, tGenome, hal, output_maf)
	command = subprocess.check_output(cmd, shell = True).decode()
	# Parse multiple sequence alignment
	for multiple_alignment in AlignIO.parse(output_maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(tGenome):
				tgtChrom, tgtStart, tgtLen, tgtStrand, tgtSize = seqrec.id.split('.', 1)[1], seqrec.annotations['start'], seqrec.annotations["size"], seqrec.annotations["strand"], seqrec.annotations["srcSize"]
				# Change coordinates following the strand
				strand = ['+' if tgtStrand == 1  else '-']
				if strand[0] == '+':
					start = tgtStart
					end = tgtStart + tgtLen
				else:
					end = tgtSize - tgtStart
					start = tgtSize - (tgtStart + tgtLen)
				if tgtChrom == tSeq:
					mypos[tgtChrom, tStart, tEnd].append([start, end])

	with open(output_qc, 'w') as quality_check:
		for (chromosome_t, start_t, end_t) in mypos:
			# Let's take 100 base pairs upstream and downstream the BUSCO coordinates to allow some flexiblity
			start_gene = int(start_t) - 100
			end_gene = int(end_t) + 100
			n = 0
			for (start, end) in mypos[chromosome_t, start_t, end_t]:
				if start in range(start_gene, end_gene) and end in range(start_gene, end_gene):
					len = end - start
					n += len
			quality_check.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(qGenome, qSeq, qStart, qEnd, tGenome, tSeq, start_t, end_t, n))
	os.remove(output_maf)





if __name__ == "__main__":
	args = parser.parse_args()
	# Generate all pairwise comparisons
	pairwise_combinations = pairwise_comparisons(args.refGenome, args.t, args.o)
	species_tol_id = get_species_name_tol_id(args.f)
	for (target, query) in pairwise_combinations:
		single_copy_orthogroups = get_single_copy_orthogroups(args.ort, query, target, main_path, species_tol_id, args.hal, args.o)


