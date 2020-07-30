#!/usr/bin/env python3

import os
import sys
import re
import shutil

#additional modules

import pyfaidx



def modifyFasta(REF,region,repeat,maxsize,contig):

	'''
	Generate synthetic chromosomes
	'''

	fastafile=pyfaidx.Fasta(REF)
	CHROM,START,END=region.split(':')[0],int(region.split(':')[1].split('-')[0]),int(region.split(':')[1].split('-')[1])
	refsequence=fastafile[CHROM][:len(fastafile[CHROM])].seq[START-1:END]

	r=re.compile(r'('  + repeat + r')\1+')
	seen=set()

	for match in r.finditer(refsequence):

		seen.add(match.group(0))

	if len(seen) == 0:

		print('[Error] Could not find  repeat in reference')
		sys.exit(1)

	else:

		min_=sum(el.count(repeat) for el in seen)
		print('[Message] ' + str(min_) + ' ' + str(repeat) + ' in reference sequence')

		number=(maxsize-min_)/contig
		seqdict=dict()

		for reps in range(min_+round(number), maxsize, round(number)):

			header='>treadmill_reef_' + region + '_' + repeat + '_' + str(reps)
			sequence=refsequence+repeat*(reps-min_)
			seqdict[header]=sequence

		fasta=''

		for key,value in seqdict.items():

			fasta+=key+'\n'+value+'\n'

		return fasta


def run(parser,args):

	'''
	Execute the code and store to FASTA
	'''

	REFIN=os.path.abspath(args.fastafile)

	try:

		with open(REFIN) as genome:

			assert(genome.readline().startswith('>'))

	except:

		print('[Error] Invalid reference FASTA file')
		sys.exit(1)

	REFOUT=os.path.abspath(args.output)

	if not os.access(os.path.dirname(REFOUT),os.W_OK):

		print('[Error] Missing write permissions on the output folder')
		sys.exit(1)

	if os.path.exists(os.path.abspath(REFOUT)):

		print('[Error] Output file already exists.')
		sys.exit(1)

	r=re.compile('.*:.*-.*')

	if r.match(args.region) is None:

		print('[Error] Wrong region format')
		sys.exit(1)

	shutil.copy2(REFIN,REFOUT)
	sequence=modifyFasta(REF,args.region,args.repeat,args.maxsize,args.contig)

	with open(REFOUT, 'a') as fout:

		fout.write(sequence)

	sys.exit(0)