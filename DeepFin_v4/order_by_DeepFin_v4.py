#!/usr/bin/env python


# Author : @cb46


import csv
from sys import argv



deepfin = argv[1]
species_list = argv[2]
output = open('DeepFin_v4_classification.tsv', 'w')



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
			if line[3] == order and line[4] == family:
				assembly_version = line[0].split('.', 1)[1]
				output.write('{}\t{}\t{}\t{}\n'.format(line[1], keys, line[5], assembly_version))

        
        
   
