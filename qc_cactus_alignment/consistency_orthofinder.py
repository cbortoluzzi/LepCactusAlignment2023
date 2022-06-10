#!/usr/bin/env python


# Author : @cb46


import os
import argparse
import subprocess
from ete3 import Tree
from pathlib import Path
from Bio import AlignIO
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Check consistency of Orthogroups genomic coordinates with cactus alignment')
parser.add_argument('--refGenome', help = 'Name of reference specie to use as query')
parser.add_argument('--list_orthogroups', help = 'List of single-copy orthogroups, one per line')
parser.add_argument('--species_list', help = 'Tab delimited species list file')
parser.add_argument('--tree', help = 'Phylogenetic tree used as guide tree in cactus')
parser.add_argument('--hal', help = 'Input hal file')
parser.add_argument('--o', help = 'Output directory')



def pairwise_comparisons(refGenome, tree, path):
	list_nodes = []
	phylo = Tree(tree, format = 1)
	# Re-root the guide tree to the outgroup
	phylo.set_outgroup('Hydropsyche_tenuis')
	for node in phylo.traverse('postorder'):
		if node.is_leaf():
			if node.name != refGenome and node.name != "Hydropsyche_tenuis":
				# Phylogenetic distance
				output_directory = Path(path, refGenome + '_vs_' + node.name)
				output_directory.mkdir(parents=True, exist_ok=True)
				list_nodes.append(node.name)
	mypair = list(map(lambda e: (e, refGenome), list_nodes))
	return mypair



def get_species_name_tol_id(species_list):
	species_d = {}
	with open(species_list) as f:
		for line in f:
			assembly, tol_id, p_class, species_name, superfamily = line.strip().split()
			species_d[species_name] = [tol_id, assembly, superfamily]
	return species_d



def get_genomic_coordinates(species):
	mygene = {}
	path = '/lustre/scratch123/tol/teams/durbin/users/cb46/whole_genome_alignment/lepidoptera/89way/version_202201/quality_check/orthofinder/proteomes/primary_transcripts_coordinates'
	pep_f = Path(path, species + '.bed')
	if Path(pep_f).is_file():
		with open(pep_f) as f:
			for line in f:
				gene_name, chromosome, start, end = line.strip().split()
				mygene[gene_name] = [chromosome, start, end]
	return mygene



def get_single_copy_orthogroups(orthogroups, query, target, gene_coordinates_query, gene_coordinates_target, species_d, hal, path):
	with open(orthogroups) as f:
		for line in f:
			orthogroup_id = line.strip()
			gene_names = get_gene_name(orthogroup_id, query, target, gene_coordinates_query, gene_coordinates_target)
			for key, value in gene_names.items():
				output_maf = Path(path, key[1] + '_vs_' + key[2], key[0] + '.maf')
				output_qc = Path(path, key[1] + '_vs_' + key[2], key[0] + '.qc')
				genome_query = key[1].lower() + '_' + species_d[key[1]][1].replace('.', 'v').replace('GCA_', 'gca')
				genome_target = key[2].lower() + '_' + species_d[key[2]][1].replace('.', 'v').replace('GCA_', 'gca')
				# Query species
				refSequence = value[0][0]
				refStart = value[0][1]
				refEnd = value[0][2]
				length = int(refEnd) - int(refStart)
				# Target species
				targetSequence = value[1][0]
				targetStart = value[1][1]
				targetEnd = value[1][2]
				maf = pairwise_alignment(refSequence, genome_query, refStart, refEnd, length, genome_target, hal, output_maf, output_qc, targetSequence, targetStart, targetEnd)

	return gene_names



def get_gene_name(orthogroup_id, query, target, gene_coordinates_query, gene_coordinates_target):
	dictionary = {}
	path = '/lustre/scratch123/tol/teams/durbin/users/cb46/whole_genome_alignment/lepidoptera/89way/version_202201/quality_check/orthofinder'
	pep_query_vs_pep_target = Path(path, 'proteomes/primary_transcripts/OrthoFinder/Results_Feb03/Orthologues/', 'Orthologues_' + query + '.pep', query + '.pep__v__' + target + '.pep.tsv')
	if Path(pep_query_vs_pep_target).is_file():
		with open(pep_query_vs_pep_target) as f:
			next(f)
			for line in f:
				line = line.strip().split()
				if orthogroup_id == line[0]:
					gene_query = line[1]
					gene_target = line[2]
					dictionary[orthogroup_id, query, target] = [gene_coordinates_query[gene_query], gene_coordinates_target[gene_target]]
	return dictionary



def pairwise_alignment(refSequence, genome_query, refStart, refEnd, length, genome_target, hal, output_maf, output_qc, targetSequence, targetStart, targetEnd):
	mypos = defaultdict(list)
	cmd = 'hal2maf --refSequence %s --refGenome %s --start %s --length %s --noAncestors --onlyOrthologs --targetGenomes %s %s %s' %(refSequence, genome_query, refStart, length, genome_target, hal, output_maf)
	try:
		command = subprocess.check_output(cmd, shell = True).decode()
	except:
		raise ValueError('hal2maf failed!')
	# Parse multiple sequence alignment
	for multiple_alignment in AlignIO.parse(output_maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(genome_target):
				tgtChrom, tgtStart, tgtLen, tgtStrand, tgtSize = seqrec.id.split('.', 1)[1], seqrec.annotations['start'], seqrec.annotations["size"], seqrec.annotations["strand"], seqrec.annotations["srcSize"]
				# Change coordinates following the strand
				strand = ['+' if tgtStrand == 1  else '-']
				if strand[0] == '+':
					start = tgtStart
					end = tgtStart + tgtLen
				else:
					end = tgtSize - tgtStart
					start = tgtSize - (tgtStart + tgtLen)
				if tgtChrom == targetSequence:
					mypos[tgtChrom, targetStart, targetEnd].append([start, end])

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
			quality_check.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(genome_query, refSequence, refStart, refEnd, genome_target, targetSequence, start_t, end_t, n))
	os.remove(output_maf)





if __name__ == "__main__":
	args = parser.parse_args()
	# Generate all pairwise comparisons
	pairwise_combinations = pairwise_comparisons(args.refGenome, args.tree, args.o)
	species_tol_id = get_species_name_tol_id(args.species_list)
	for (target, query) in pairwise_combinations:
		# Loop through each comparison
		genomic_coordinates_query = get_genomic_coordinates(query)
		genomic_coordinates_target = get_genomic_coordinates(target)
		single_copy_orthogroups = get_single_copy_orthogroups(args.list_orthogroups, query, target, genomic_coordinates_query, genomic_coordinates_target, species_tol_id, args.hal, args.o)

	
