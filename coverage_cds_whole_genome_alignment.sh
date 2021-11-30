#!/bin/bash


# Author: @cb46


# Export path
export PATH=/software/team118/maf_stream/target/release:$PATH



if [[ -z $1 ]] ; then
  echo "Usage: ./coverage_cds_whole_genome_alignment.sh <GFF> <HAL> <species.name> "
  exit -1
fi



GFF=$1
HAL=$2
SPECIES=$3
DATE="$(date +'%Y%m%d')"
DIR="coverage_cds_"$DATE


mkdir -p $DIR


echo "./coverage_cds_whole_genome_alignment.sh $GFF $HAL $SPECIES"


# Extract from annotation the genomic coordinates of all coding sequences
cat $GFF | grep -v '#' | awk 'OFS="\t"{if($3=="CDS"){print $1,$4,$5}}' | sortBed | mergeBed -i stdin | awk 'OFS="\t"{print $1,$2,$3,$3-$2}' > $SPECIES.cds.bed


# Get alignment depth at each position in the multiple sequence alignment
cat $SPECIES.cds.bed | while read contig start end length
do
  echo "hal2maf --noAncestors --refSequence $contig --start $start --length $length --refGenome $SPECIES --append $HAL $SPECIES.maf"
  hal2maf --noAncestors --refSequence $contig --start $start --length $length --refGenome $SPECIES --append $HAL $SPECIES.maf
done


echo "maf_stream coverage $SPECIES $SPECIES.maf $SPECIES.coverage.cds.txt"
maf_stream coverage $SPECIES $SPECIES.maf $SPECIES.coverage.cds.txt


plot_coverage.py $SPECIES.coverage.cds.txt



