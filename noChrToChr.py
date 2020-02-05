#!/usr/bin/env python


#Hasiba Asma 
#21st Aug 2018


import os 
import sys

file=sys.argv[1]

with open(file,'r') as infile,open('with_chr'+file,'w') as outfile:
	for line in infile:
		line="chr"+line
		outfile.write(line)
		#col=line.split('\t')
		#chrName=col[0]

		#line=line.strip('chr')
