#!/usr/bin/env python


# Author : @cb46


import argparse
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Define intronic regions from coding sequence coordinates')
parser.add_argument('--cds', help = 'Tab delimited file with coding sequence coordinates')


def define_intronic_regions(bed_f, output_d, refGenome):
        cds = defaultdict(list)
        with open(bed_f) as f:
                for line in f:
                        chromosome, start, end, feature, strand, gene = line.strip().split()
                        start = int(start)
                        end = int(end)
                        cds[chromosome, gene, strand].append([start, end])

        output_f = Path(output_d, refGenome + '.gff3.intron.bed')
        with open(output_f, 'w') as out:
                for gene in cds:
                        # We will consider only genes with more than one coding sequence
                        if len(cds[gene]) > 1:
                                list_s = [i[0] for i in cds[gene]]
                                list_e = [i[1] for i in cds[gene]]
                                intron = list(zip(list_e, list_s[1:]))
                                for genomic_coordinates in intron:
                                        out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene[0], genomic_coordinates[0], genomic_coordinates[1], 'intron', gene[2], gene[1]))




if __name__ == "__main__":
        args = parser.parse_args()
        # Create output directory if it doesn't exist
        path = Path(args.cds).parent
        path.mkdir(parents=True, exist_ok=True)
        refGenome = Path(args.cds).stem.split('.')[0]
        intron = define_intronic_regions(args.cds, path, refGenome)

