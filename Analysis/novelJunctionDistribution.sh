#!/bin/bash

db='SpliceJunction.db'
home='/home/dennis/Github/MendelianRNA-seq-DB/Analysis'
list=`cat control.list`
count=`cat count.list`

python3 $home/FilterSpliceJunctions.py --sample 10-1-M.bam 5

for sample in $list; do

	python3 $home/FilterSpliceJunctions.py --delete $sample.bam
	python3 $home/FilterSpliceJunctions.py --sample 10-1-M.bam 5

done

for file in $count; do

	wc -l $file
done