#!/usr/bin/env python3

import os
import sys
import re

#additional modules

import pyfaidx



def run(parser,args):

	'''
	Execute the code and store to FASTA
	'''

	REF=os.path.abspath(args.fastafile)

	try:

		with open(REF) as genome:

			assert(genome.readline().startswith('>'))

	except:

		print('[Error] Invalid reference FASTA file')
		sys.exit(1)

	FOUT=os.path.abspath(args.output)

	if not os.access(os.path.dirname(FOUT),os.W_OK):

		print('[Error] Missing write permissions on the output folder')
		sys.exit(1)


	r=re.compile('.*:.*-.*')

	if r.match(args.region) is None:

		print('[Error] Wrong region format')
		sys.exit(1)

	sys.exit(0)