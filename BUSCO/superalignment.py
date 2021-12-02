#!/usr/bin/python


# Author: @cb46


import argparse
from pathlib import Path


parser = argparse.ArgumentParser(description = 'Generate supermatrix for RAxML')
parser.add_argument('-i', help = 'Name input directory')
parser.add_argument('-o', help = 'Name output directory')



def generate_supermatrix(in_folder, out_folder):
	data = {}
	key = None
	files = list(Path(in_folder).rglob('*.trimal'))
	for f in files:
		for line in open(f, 'r'):
			if line.startswith('>'):
				key = line.strip()[1:]
				if key in data:
					pass
				else:
					data.update({key: ''})
			else:
				data[key] += line.strip()

	output_file = 'supermatrix.aln.faa'
	output = Path(out_folder, output_file)
	with open(output, 'w') as out:
		for key in data:
			out.write('>%s\n' % key)
			out.write('%s\n' % data[key])



if __name__ == '__main__':
	args = parser.parse_args()
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	supermatrix = generate_supermatrix(args.i, p)

	
	
	
