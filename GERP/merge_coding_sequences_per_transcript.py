#!/usr/bin/env python


# Author : @cb46


import argparse
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Obtain unique non-overlapping instances of each feature type')
parser.add_argument('--gff3', help = 'Annotation in GFF3 format')
parser.add_argument('--feature', help = 'Feature type')
parser.add_argument('--refGenome', help = 'Reference genome')


def parse_annotation(gff_f, feature):
	mygff = defaultdict(list)
	with open(gff_f) as f:
		for line in f:
			if not line.startswith('#'):
				contig, ensembl, feature_type, start, end, dot, strand, value, info = line.strip().split()
				start, end = int(start), int(end)
				if feature_type == "gene":
					gene_id, gene_biotype, gene_id, gene_version = info.split(';')
					gene = gene_id.split('=')[1]
					gene_contig, gene_start, gene_end, gene_strand = contig, start, end, strand
					mygff[gene_contig, "gene", gene_start, gene_end, gene_strand, gene]
				elif feature_type == feature:
					feature_start, feature_end = start, end
					mygff[gene_contig, "gene", gene_start, gene_end, gene_strand, gene].append([feature_start, feature_end])
	return mygff



def unique_nonoverlapping_instances(mygff, feature, refGenome, path):
	for gene in mygff:
		output_file = Path(path, refGenome+'.'+gene[0]+'.CDS.gff3')
		with open(output_file, 'a') as out:
			intervals = mygff[gene]
			if intervals:
				intervals.sort(key=lambda interval: interval[0])
				merged = [intervals[0]]
				for current in intervals:
					previous = merged[-1]
					if current[0] <= previous[1]:
						previous[1] = max(previous[1], current[1])
					else:
						merged.append(current)
				for item in merged:
					contig = refGenome+'.'+gene[0]
					info = "Parent=gene"+gene[5]
					out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(contig, 'ensembl', feature, item[0], item[1], '.', gene[4], '.', info))




if __name__ == "__main__":
	args = parser.parse_args()
	path = Path(args.gff3).parent
	select_features = parse_annotation(args.gff3, args.feature)
	unique_elements = unique_nonoverlapping_instances(select_features, args.feature, args.refGenome, path)


