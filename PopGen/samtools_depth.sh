#!/bin/bash



# Author: @cb46


for bam in bam/*.bam
do
        echo "samtools depth $bam"
        /software/team118/samtools/1.11/bin/samtools depth $bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $bam".cov"
done

