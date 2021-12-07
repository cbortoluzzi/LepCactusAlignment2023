#!/usr/bin/env python


# Author : @cb46


import csv
from sys import argv



def deepfin_classification_bony_fishes(deepfin, species_list, output, missclassified):
	mydict = {}
	with open(deepfin) as csvfile:
		f = csv.reader(csvfile, delimiter=',')
		headers = next(f, None)
		header = '\t'.join(headers[1:-4])
		output.write('{}\t{}\t{}\t{}\n'.format('ToL', header, 'Species', 'Assembly'))
		next(f)
		for line in f:
			key = "\t".join(line[1:-4])
			mydict[key] = []

	for keys in mydict:
		with open(species_list) as f:
			next(f)
			for line in f:
				line = line.strip().split('\t')
				key = keys.strip().split('\t')
				order, family = key[-4], key[-1]
				assembly_version = line[0].split('.', 1)[1]
				# Let's first check that the order and the family name in the DeepFin classification are equal to the ones reported in our tsv file
				if line[3] == order and line[4] == family:
					output.write('{}\t{}\t{}\t{}\n'.format(line[1], keys, line[5], assembly_version))
				# If that's not the case, check that the famly name in the DeepFin classification is equal to that reported in the tsv file
				elif line[4] == family:
					output.write('{}\t{}\t{}\t{}\n'.format(line[1], keys, line[5], assembly_version))
					missclassified.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(line[5], line[3], line[4], "-->" , order, family))
				
				
				
if __name__ == "__main__":
	deepfin, species_list = argv[1], argv[2]
	output = open('DeepFin_v4_classification.tsv', 'w')
	missclassified = open('misclassified_species.tsv', 'w')
	deepfin_classification = deepfin_classification_bony_fishes(deepfin, species_list, output, missclassified)
	output.close()
	missclassified.close()
	
	
