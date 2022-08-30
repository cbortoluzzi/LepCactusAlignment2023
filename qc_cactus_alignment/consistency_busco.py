#!/usr/bin/env python


# Author : @cb46


import os
import argparse
import itertools
import subprocess
from ete3 import Tree
from pathlib import Path
from Bio import AlignIO
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Check consistency of single copy BUSCO genes genomic coordinates with cactus alignment')
parser.add_argument('--refGenome', help = 'Name of species to use as query')
parser.add_argument('--list_genes', help = 'List of single-copy BUSCO genes, one per line')
parser.add_argument('--species_list', help = 'A tab delimited species list file')
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
			# Remove the outgroup from the pairwise comparisons
			if node.name != refGenome and node.name != "Hydropsyche_tenuis":
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
			species_d[species_name] = [tol_id, assembly]
	return species_d



def busco_genes(list_genes):
	genes = []
	with open(list_genes) as f:
		for line in f:
			busco = line.strip()
			genes.append(busco)
	return genes



def change_assembly_coordinates(query, target, species_d):
	query_d, target_d = {}, {}
	main_path = '/lustre/scratch123/tol/projects/lepidoptera/freeze/2021-06-21/data/insects'
	alternative_path = '/lustre/scratch123/tol/teams/durbin/users/cb46/darwin/lepidoptera/busco.v5'
	# We need to change the chromosome nomenclature to match the one in the cactus alignment
	# For example, HG992306.1 is chromosome '1' in Tinea trinotella
	change_q = change_coordinates(main_path, alternative_path, query, species_d[query][0], species_d[query][1], query_d)
	change_t = change_coordinates(main_path, alternative_path, target, species_d[target][0], species_d[target][1], target_d)
	return change_q, change_t



def change_coordinates(main_path, alternative_path, species_name, species_tol_id, species_assembly, dictionary):
	assembly_report = Path(main_path, species_name, 'assembly', 'release', species_tol_id, 'insdc', species_assembly+'_assembly_report.txt')
	# Check that assembly report file exists
	if not Path(assembly_report).is_file():
		assembly_report = Path(alternative_path, species_name, species_assembly+'_assembly_report.txt')
	with open(assembly_report) as f:
		for line in f:
			if not line.startswith('#'):
				line = line.strip().split()
				if line[1] == 'assembled-molecule':
					dictionary[line[4]] = line[2]
				elif line[1] == 'unplaced-scaffold':
					dictionary[line[4]] = line[0]
				elif line[1] == 'unlocalized-scaffold':
					dictionary[line[4]] = line[0]
	return dictionary



def get_busco_coordinates(query, target, species_d, genes, change_q, change_t, hal, path):
	main_path = '/lustre/scratch123/tol/projects/lepidoptera/freeze/2021-06-21/data/insects'
	alternative_path = '/lustre/scratch123/tol/teams/durbin/users/cb46/darwin/lepidoptera/busco.v5'
	for gene in genes:
		coordinates_q = parse_fasta(main_path, alternative_path, query, species_d[query][0], gene)
		coordinates_t = parse_fasta(main_path, alternative_path, target, species_d[target][0], gene)
		with open(coordinates_q) as gene_q, open(coordinates_t) as gene_t:
			for line_q, line_t in zip(gene_q, gene_t):
				if line_q.startswith('>') and line_t.startswith('>'):
					chrom_q = line_q.strip().split(':')[0].replace('>', '')
					chrom_t = line_t.strip().split(':')[0].replace('>', '')
					start_end_q = line_q.strip().split(':')[1]
					start_end_t = line_t.strip().split(':')[1]
					start_q, end_q = start_end_q.split('-')
					start_t, end_t = start_end_t.split('-')
					# Change chromosome same so that the one in the alignment corresponds to that from BUSCO
					chromosome_q = change_q[chrom_q]
					chromosome_t = change_t[chrom_t]
					genome_query = query.lower()+'_'+species_d[query][1].lower().replace('_', '').replace('.', 'v')
					genome_target = target.lower()+'_'+species_d[target][1].lower().replace('_', '').replace('.', 'v')
					# Obtain multiple sequence alignment for each gene of query and target species
					maf = multiple_alignment_format(genome_query, chromosome_q, start_q, end_q, genome_target, chromosome_t, start_t, end_t, gene, hal, path)



def parse_fasta(main_path, alternative_path, species_name, species_tol_id, gene):
	# Let's get the genomic coordinates (i.e. chromosomem, start, end) of each single-copy BUSCO gene
	faa = Path(main_path, species_name, 'analysis', species_tol_id, 'busco', 'lepidoptera_odb10_metaeuk', 'run_lepidoptera_odb10', 'busco_sequences', 'single_copy_busco_sequences', gene+'.faa')
	# Check if file exists
	if not Path(faa).is_file():
		faa = Path(alternative_path, species_name, 'lepidoptera_odb10_metaeuk', 'run_lepidoptera_odb10', 'busco_sequences', 'single_copy_busco_sequences', gene+'.faa')
	return faa



def multiple_alignment_format(genome_query, chromosome_q, start_q, end_q, genome_target, chromosome_t, start_t, end_t, gene, hal, path):
	mypos = defaultdict(list)
	len_q = int(end_q) - int(start_q)
	output_maf = Path(path, query + '_vs_' + target, gene + '.maf')
	output_qc = Path(path, query + '_vs_' + target, gene + '.qc')
	# Run hal2maf to obtain an alignment in multiple alignment format (MAF)
	cmd = 'hal2maf --refSequence %s --refGenome %s --start %s --length %s --noAncestors --onlyOrthologs --targetGenomes %s %s %s' %(chromosome_q, genome_query, start_q, len_q, genome_target, hal, output_maf)
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
				if tgtChrom == chromosome_t:
					mypos[tgtChrom, gene, start_t, end_t].append([start, end])

	with open(output_qc, 'w') as quality_check:
		for (chromosome, gene, start_g, end_g) in mypos:
			# Let's take 100 base pairs upstream and downstream the BUSCO coordinates to allow some flexiblity
			start_gene = int(start_g) - 100
			end_gene = int(end_g) + 100
			n = 0
			for (start, end) in mypos[chromosome, gene, start_g, end_g]:
				if start in range(start_gene, end_gene) and end in range(start_gene, end_gene):
					len = end - start
					n += len
			quality_check.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(genome_query, chromosome_q, start_q, end_q, genome_target, chromosome_t, start_g, end_g, n))
	# We can now remove the MAF file
	os.remove(output_maf)




if __name__ == "__main__":
	args = parser.parse_args()
	# Generate all possible pairwise comparisons between the query and the target species
	pairwise_combinations = pairwise_comparisons(args.refGenome, args.tree, args.o)
	species_tol_id = get_species_name_tol_id(args.species_list)
	list_busco_orthogroups = busco_genes(args.list_genes)
	# Iterate over each pairwise comparison and each single-copy BUSCO gene
	for (target, query) in pairwise_combinations:
		# Update genomic coordinates of BUSCO genes using the assembly report file
		assembly_coordinates_query, assembly_coordinates_target = change_assembly_coordinates(query, target, species_tol_id)
		# Obtain an alignment in multiple alignment format (MAF) for each BUSCO gene and check consistency
		evaluate_consistency = get_busco_coordinates(query, target, species_tol_id, list_busco_orthogroups, assembly_coordinates_query, assembly_coordinates_target, args.hal, args.o)

		
