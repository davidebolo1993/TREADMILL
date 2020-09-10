#!/usr/bin/env python3

#python 3 standard library

import sys
import os
import pickle
import subprocess


#additional modules


CS_CPP=os.path.abspath(os.path.dirname(__file__) + '/consensus')

def ParseGroups(BIN,OUT,match,mismatch,gapopen,gapextend):

	binin=open(BIN,'rb')
	dictR = pickle.load(binin)
	binin.close()

	for keyR in dictR.keys():

		OUTR=os.path.abspath(OUT + '/' + keyR)

		if not os.path.exists(OUTR):

			os.makedirs(OUTR)

		refsequence=list(dictR[keyR]['group0']['reference'].keys())[0]

		with open(os.path.abspath(OUTR + '/r.tmp.fa'), 'w') as fr:

			fr.write('>reference\n' + refsequence + '\n')

		for i,keyG in enumerate(dictR[keyR].keys()):

			if keyG != 'group0':

				counter=0 #do not exceed an arbitraty treshold, otherwise consensus sequence computation becomes slow


				with open(os.path.abspath(OUTR + '/a'+str(i+1)+'.tmp.fa'), 'a') as fa:

					for keyS in dictR[keyR][keyG].keys():

						counter +=1

						if counter < 100: #arbitrary treshold

							fa.write('>'+keyS+'\n'+list(dictR[keyR][keyG][keyS].keys())[0]+'\n')

						else:

							continue

				with open(os.path.abspath(OUTR + '/a'+str(i+1)+'.cs.fa'), 'w') as cs:

					subprocess.call([CS_CPP, str(match), str(mismatch), str(gapopen), str(gapextend), os.path.abspath(OUTR + '/r.tmp.fa'), os.path.abspath(OUTR + '/a'+str(i+1)+'.tmp.fa')], stdout=cs, stderr=open(os.devnull, 'wb'))

				os.remove(os.path.abspath(OUTR + '/a'+str(i+1)+'.tmp.fa')) #clean-up


def run(parser,args):

	'''
	Execute the code and write variants to BCF
	'''

	BIN=os.path.abspath(args.input)

	if not os.path.isfile(BIN):

		print('[Error] Invalid BIN file')
		sys.exit(1)

	OUT=os.path.abspath(args.output)

	match=args.match
	mismatch=args.mismatch
	gapopen=args.gapopen
	gapextend=args.gapextend

	if not os.path.exists(OUT):

		os.makedirs(OUT)

	ParseGroups(BIN,OUT,match,mismatch,gapopen,gapextend)