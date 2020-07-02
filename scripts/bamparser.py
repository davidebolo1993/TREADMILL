#!/usr/bin/env python3

#python 3 standard library

import sys
import os
import re
import json
from collections import defaultdict

#additional modules

import pysam

#global

S_dict=dict()

S_dict['BAM_CSOFT_CLIP']=0
S_dict['BAM_CMATCH']=0
S_dict['BAM_CINS']=0
S_dict['BAM_CDEL']=0
S_dict['BAM_CDIFF']=0
S_dict['BAM_UNMAPPED']=0
S_dict['BAM_SUPP']=0
S_dict['BAM_SEC']=0



def tuptodict(cigartuples):

	'''
	Convert tuples from read.cigartuples to dict and sum values from same keys
	'''

	Cdict = defaultdict(int)

	for key, value in cigartuples:

		Cdict[key] = Cdict[key]+value

	return Cdict


def parse(BAM,query): 

	'''
	Parse BAM to get statistics about mapped and unmapped reads, optionally focusing on a specific region/chromosome
	'''

	global S_dict

	bamfile=pysam.AlignmentFile(BAM, "rb")

	if not query:

		it=bamfile.fetch()

	else:

		if len(query.split(":")) == 1:

			it=bamfile.fetch(query.split(":")[0])

		else:

			it=bamfile.fetch(query.split(":")[0],int(query.split(":")[1].split("-")[0]),int(query.split(":")[1].split("-")[1]))

	for read in it:

		if read.has_tag('MD'):

			if not read.is_unmapped:

				if not read.is_supplementary:

					if not read.is_secondary:

						Cdict=tuptodict(read.cigartuples)
						MD=read.get_aligned_pairs(with_seq=True)

						S_dict['BAM_CSOFT_CLIP']+=Cdict[4]
						S_dict['BAM_CINS']+=Cdict[1]
						S_dict['BAM_CDEL']+=Cdict[2]
						S_dict['BAM_CMATCH']+=sum(1 for x,y,z in MD if x is not None and z is not None and z[0].isupper())
						S_dict['BAM_CDIFF']+=sum(1 for x,y,z in MD if x is not None and z is not None and z[0].islower())

					else:

						S_dict['BAM_SEC']+=read.infer_query_length()

				else:

					S_dict['BAM_SUPP']+=read.infer_query_length()

			else:

				S_dict['BAM_UNMAPPED']+=read.infer_query_length()

		else:

			print('[Error] Missing required MD tag')
			bamfile.close()
			sys.exit(1)

	bamfile.close()

	return S_dict


def main():

	'''
	Execute the code and dump json to stdout
	'''

	if len(sys.argv) > 2 or len(sys.argv) == 1:

		print("[Usage] python bamparser.py <input.bam> <chr:start-end>")
		sys.exit(1)

	elif len(sys.argv) == 2:

		BAM=os.path.abspath(sys.argv[1])
		query=False

	else:

		BAM=os.path.abspath(sys.argv[1])
		query=sys.argv[2]

	if os.path.isfile(BAM):

		S_dict=parse(BAM,query)
		json.dump(S_dict, sys.stdout)

	else:

		print('[Error] Invalid BAM file')
		sys.exit(1)


if __name__=='__main__':

	main()
