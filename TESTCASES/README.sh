#!/bin/bash

#Run this script to see an example execution of the BLAST-QC program.
#The output of each execution is placed into a separate subdirectory of the TESTCASES directory.

#Most basic call- filename and type of BLAST (in this case blastx) run are required parameters.
#will result in 3 files: one header with the <Query_def> of the matching hits, one with info on the hits that did not 
#match input (in this case none) and one with info on hits that match the input arguments. 
python ../BLAST-QC.py --input sampleResults.xml --type n --output one.out/example_one.out --order e

#replicating a `-max_target_seqs` parameter (without all the confusion) note the low level of detail in the hit's definition.
python ../BLAST-QC.py -i sampleResults.xml -t n -o two.out/example_two.out --number 1 -or e

#There! using a range value makes finding a good hit with acceptable definition info easier.
python ../BLAST-QC.py -i sampleResults.xml -t n -o three.out/example_three.out -n 2 --erange .0005 -or e

#Thresholds for lots of important values.
python ../BLAST-QC.py -i sampleResults.xml -t n -o four.out/example_four.out -n 1 -or e --evalue .03 --bitscore 35 --description 1 --identity 60

