#!/usr/bin/env python


#Hasiba Asma 
#21st Aug 2018


import os 
import sys

file=sys.argv[1]

with open(file,'r') as infile,open('without_chr'+file,'w') as outfile:
	for line in infile:
		line=line.strip('chr')
		outfile.write(line)
		#col=line.split('\t')
		#chrName=col[0].strip('chr')
