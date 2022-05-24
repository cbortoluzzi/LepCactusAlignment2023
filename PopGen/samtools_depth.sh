#!/bin/bash



# Author: @cb46


for BAM in bam/*.bam
do
        echo "samtools depth $BAM"
        /software/team118/samtools/1.11/bin/samtools depth $BAM | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $BAM".cov"
done

