#!/bin/bash

FILE=$1

awk -v OFS='\t' '{print "chr"$1, $2, $3}' $FILE > "gb_"$FILE
