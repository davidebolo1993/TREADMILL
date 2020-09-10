#!/usr/bin/env python3

#python 3 standard library

import sys
import os
import pickle
import subprocess
from collections import OrderedDict

#additional modules

import editdistance

CS_CPP=os.path.abspath(os.path.dirname(__file__) + '/consensus')


def GTLH(alleles,coverage,error):

	'''
	
	'''





def ParseGroups(BIN,OUT,match,mismatch,gapopen,gapextend,treshold):

	'''
	'''

	binin=open(BIN,'rb')
	dictR = pickle.load(binin)
	binin.close()

	for keyR in dictR.keys():

		OUTR=os.path.abspath(OUT + '/' + keyR)

		if not os.path.exists(OUTR):

			os.makedirs(OUTR)

		refsequence=dictR[keyR]['reference']

		with open(os.path.abspath(OUTR + '/r.tmp.fa'), 'w') as fr:

			fr.write('>reference\n' + refsequence + '\n')

		#assuming diploidy

		PHom=1.0 - (2*dictR[keyR]['error'])
		PHet= .5-dictR[keyR]['error']
		listA=list()

		for i,keyG in enumerate(dictR[keyR].keys()):

			if keyG != 'reference' and keyG != 'error' and keyG != 'coverage':

				dictR[keyR][keyG] = OrderedDict(sorted(dictR[keyR][keyG].items(), key=lambda x:x[1][1])) #sort by quality. If the number of sequences is > 100, this allows to exclude low-quality sequences from consensus computation
				counter=0

				with open(os.path.abspath(OUTR + '/a'+str(i+1)+'.tmp.fa'), 'w') as fa:
				
					for keyS in dictR[keyR][keyG].keys():

						counter+=1

						if counter <= 100: #arbitrary treshold

							fa.write('>'+keyS+'\n'+dictR[keyR][keyG][keyS][0]+'\n')

						else:

							continue				

				with open(os.path.abspath(OUTR + '/a'+str(i+1)+'.cs.fa'), 'w') as cs:

					subprocess.call([CS_CPP, str(match), str(mismatch), str(gapopen), str(gapextend), os.path.abspath(OUTR + '/r.tmp.fa'), os.path.abspath(OUTR + '/a'+str(i+1)+'.tmp.fa')], stdout=cs)

				os.remove(os.path.abspath(OUTR + '/a'+str(i+1)+'.tmp.fa')) #clean-up

				with open(os.path.abspath(OUTR + '/a'+str(i+1)+'.cs.fa')) as csin:

					for line in csin:

						if line[0] != '>':

							listA.append((line.rstrip(), 1-editdistance.eval(line.rstrip(),refsequence)/max(len(line.rstrip()),len(refsequence)), counter))

				os.remove(os.path.abspath(OUTR + '/a'+str(i+1)+'.cs.fa')) #clean-up

		sortlistA=sorted(listA, key=lambda x:x[1], reverse=True)






def run(parser,args):

	'''
	Execute the code and write variants to BCF
	'''

	BIN=os.path.abspath(args.input)

	if not os.path.isfile(BIN):

		print('[Error] Invalid BIN file')
		sys.exit(1)

	OUT=os.path.abspath(args.output)

	if not os.path.exists(OUT):

		os.makedirs(OUT)

	ParseGroups(BIN,OUT,args.match,args.mismatch,args.gapopen,args.gapextend,args.similarity)



