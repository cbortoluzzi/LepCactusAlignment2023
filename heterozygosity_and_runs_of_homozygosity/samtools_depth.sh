#!/bin/bash



# Author: @cb46


bam=$1


/software/team118/samtools/1.11/bin/samtools depth $bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $bam".cov"


