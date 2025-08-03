#!/bin/bash

modkit extract full -q 500 -i 5000 -t 10 --force --kmer-size $2 --mapped-only --region $5 --reference $1 $3 $4

#$1 Reference
#$2 Kmer Size
#$3 In Bam
#$4 Outpath
#$5 Chrom