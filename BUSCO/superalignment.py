#!/usr/bin/python


# Author: @cb46


from sys import argv
from pathlib import Path


folder = argv[1]
data = {}
key = None



files = list(Path(folder).rglob('*.trimal'))
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
output = Path(folder, output_file)
with open(output, 'w') as out:
	for key in data:
		out.write('>%s\n' % key)
		out.write('%s\n' % data[key])
