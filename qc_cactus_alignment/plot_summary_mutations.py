#!/usr/bin/env python


# Author: @cb46


import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(description = 'Plot mutation events in cactus alignment')
parser.add_argument('--count', help = 'Summary table of mutation events in cactus alignemt')
parser.add_argument('--species_list', help = 'A tab delimited species list file')
parser.add_argument('--o', help = 'Output directory')



def get_genome_species_name(species_list):
        myspecies = {}
        with open(species_list) as f:
                for line in f:
                        assembly_id, assembly_name, p_class, order, family, species, iucn_red_list, project, group, genome = line.strip().split('\t')
                        myspecies[genome] = species
        return myspecies



def summary_mutation_events_alignment(count_f, myspecies, path):
        mymut = {}
        with open(count_f) as f:
                next(f)
                next(f)
                for line in f:
                        line = line.strip().split()
                        if line:
                                # Skip ancestors
                                if 'Anc' not in line[0] and 'Total' not in line[0] and 'Average' not in line[0]:
                                        genomeName = line[0].replace(',', '')
                                        species_name = myspecies[genomeName]
                                        mut = [int(i.replace(',', '')) for i in line[5:]]
                                        substitutions, transitions, transversions, insertions, deletions, inversions, duplications, transpositions, other = mut[0], mut[1], mut[2], mut[8], mut[10], mut[12], mut[14], mut[16], mut[18]
                                        mymut[species_name] = [substitutions, transitions, transversions, insertions, deletions, inversions, duplications, transpositions]
        
        df = pd.DataFrame.from_dict(mymut, orient="index", columns = ['Substitutions', 'Transitions', 'Transversions', 'Insertions', 'Deletions', 'Inversions', 'Duplications', 'Transpositions'])
        log_df = df.apply(lambda x: np.log10(x) if np.issubdtype(x.dtype, np.number) else x)
        plt.rcParams.update({'font.size': 4})
        fig = plt.figure(figsize=(10, 3))
        figure = Path(path, 'summary_mutations_in_alignment.pdf')
        log_df.plot.bar(stacked = True)
        plt.ylabel('Log10 number of mutations in the alignment')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(figure, dpi = 500, bbox_inches = 'tight')



if __name__ == "__main__":
        args = parser.parse_args()
        species_name = get_genome_species_name(args.species_list)
        summary_mutations = summary_mutation_events_alignment(args.count, species_name, args.o)


