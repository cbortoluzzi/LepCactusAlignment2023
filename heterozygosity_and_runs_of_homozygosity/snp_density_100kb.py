#!/usr/bin/env python



# Author : @cb46



import vcf
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict
from statistics import mean


parser = argparse.ArgumentParser(description = 'Calculate SNP density in a 100-kb window')
parser.add_argument('--vcf', help = 'VCF file')
parser.add_argument('--bam', help = 'BAM file')
parser.add_argument('--w', help = 'Window size [default = 100000]', type = int, default = 100000)




def sequences_bam(bam_f):
	mygenome = {}
	command = 'samtools idxstats %s | cut -f 1,2' %(bam_f)
	cmd = subprocess.check_output(command, shell = True).decode()
	outcmd = cmd.split('\n')
	for line in outcmd:
		if line:
			chromosome, length = line.strip().split()
			length = int(length)
			try:
				if str(chromosome) == 'W' or str(chromosome) == 'Z' or isinstance(int(chromosome), int):
					mygenome[chromosome] = length
			except ValueError:
				continue
				#print (f"Contig {chromosome} is not considered")
	return mygenome



def calculate_snp_density(mygenome, vcf_f, window, output_file):
	mydict = defaultdict(list)
	for chromosome in mygenome:
		seq_len = mygenome[chromosome]
		for i in range(0, seq_len, window):
			density_per_100kb = density_distribution(vcf_f, chromosome, i, i + window, mydict, output_file)
	diff_list = []
	distance = []
	for i in density_per_100kb:
		list_n = density_per_100kb[i]
		for x, y in zip(list_n[0::], list_n[1::]):
			diff_list.append(y-x)
		distance.append(mean(diff_list))
	print ("Average distance between SNPs:", mean(distance))



def density_distribution(vcf_f, chromosome, start, end, mydict, output_file):
	vcf_reader = vcf.Reader(filename=vcf_f)
	nhet = 0
	with open(output_file, 'a') as out:
		for record in vcf_reader.fetch(chromosome, start, end):
			nhet += record.num_het
			for call in record.samples:
				if call['GT'] == '0/1':
					mydict[chromosome].append(record.POS)
		out.write('{}\t{}\t{}\t{}\n'.format(chromosome, start, end, nhet))
	return mydict




if __name__ == "__main__":
	args = parser.parse_args()
	output_f = Path(args.vcf).stem + '.SNPdensity100kb.bed'
	path = Path(args.vcf).parents[0]
	output_file = Path(path, output_f)
	genome = sequences_bam(args.bam)
	snp_density = calculate_snp_density(genome, args.vcf, args.w, output_file)

