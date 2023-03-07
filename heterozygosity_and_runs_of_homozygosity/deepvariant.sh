#!/bin/bash


# Author: @cb46

if [ -z $1 ]; then
  echo "Usage: ./deepvariant.sh <assembly in fasta format> <alignment in bam format> <output in vcf format>"
  exit -1
fi


fasta=$1
bam=$2
vcf=$3


# Run deepvariant
singularity exec -B /lustre deepvariant.simg /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref $fasta --reads=$bam--output_vcf=$vcf

