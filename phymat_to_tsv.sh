#!/bin/bash

#$1 input phylip distance matrix
#$2 output tab separated csv with column and row names

printf '\t' > $2; tail -n +2 $1 | cut -d' ' -f1 | tr '\n' '\t'  | sed 's/\t*$//'>> $2; printf '\n' >> $2; tail -n +2 $1 | sed 's/  /\t/g' | awk 'NF' >> $2
