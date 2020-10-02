#!/usr/bin/env python3

#python 3 standard library

import sys
import os
import re
import pickle
import datetime
import subprocess
import shutil
import gzip
from collections import OrderedDict,Counter
from itertools import combinations_with_replacement

#additional modules

import editdistance
import numpy as np
import pandas as pd
import pysam

CS_CPP=os.path.abspath(os.path.dirname(__file__) + '/consensus')

def RegExF(sequence,motif):

	'''
	Find repeats and interruptions given a known repeated motif and the sequence
	'''

	dictI=Counter()

	occurences=[(r.start(0), r.end(0)) for r in re.finditer(motif, sequence)]

	for i in range(len(occurences)-1):

		if occurences[i][1] != occurences[i+1][0]:

			dictI[sequence[occurences[i][1]:occurences[i+1][0]]]+=1

	return len(occurences), ','.join(dictI.keys())+':'+','.join(str(x) for x in dictI.values())


def VCFH(ctgs,BIN):

	'''
	Write VCF header
	'''

	#genmeral informations

	fileformat='##fileformat=VCFv4.2\n'
	filedate='##fileDate=' + ''.join(str(datetime.date.today()).split('-')) + '\n'
	source='##source=TREADMILL\n'
	contigs=''.join(x for x in ['##contig=<ID=' + x + '>\n'  for x in ctgs])

	#INFO field

	end='##INFO=<ID=END,Number=1,Type=Integer,Description="Repetition end">\n'
	motif='##INFO=<ID=MOTIF,Number=1,Type=String,Description="Repeat motif">\n'
	refrep='##INFO=<ID=RN,Number=1,Type=Integer,Description="Number of repeats in the reference sequence">\n'
	refint='##INFO=<ID=RI,Number=1,Type=String,Description="Interruptions in the reference sequence">\n'
	allelerefseq='##INFO=<ID=RALS,Number=1,Type=String,Description="Reference allele sequence (if present)">\n'
	alleledistance='##INFO=<ID=RALD,Number=1,Type=Float,Description="Sequence similarity (edit distance-based) between reference sequence and reference allele (if present)">\n'
	allelerefnum='##INFO=<ID=RALN,Number=1,Type=Integer,Description="Number of repeats in the reference allele (if present)">\n'
	allelerefint='##INFO=<ID=RALI,Number=1,Type=String,Description="Interruptions in the reference allele (if present)">\n'
	altdist='##INFO=<ID=ALTD,Number=1,Type=Float,Description="If a reference allele is not present, the highest similarity score (edit distance-based) between the alternative allele/s and the reference sequence">\n'
	allele1num='##INFO=<ID=AL1N,Number=1,Type=Integer,Description="Number of repeats in the 1st alternative allele (if present)">\n'
	allele1int='##INFO=<ID=AL1I,Number=1,Type=String,Description="Interruptions in the 1st alternative allele (if present)">\n'
	allele2num='##INFO=<ID=AL2N,Number=1,Type=Integer,Description="Number of repeats in the 2nd alternative allele (if present)">\n'
	allele2int='##INFO=<ID=AL2I,Number=1,Type=String,Description="Interruptions in the 2nd alternative allele (if present)">\n'

	#FORMAT field

	gt='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
	gl='##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled likelihoods for genotypes 0/0, 0/1, 1/1 or 0/0, 0/1, 1/1, 0/2, 1/2, 2/2>"\n'
	dp='##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
	ad='##FORMAT=<ID=AD,Number=G,Type=Integer,Description="Depth of alleles>"\n'
	head='#CHROM' + '\t' + 'POS' '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' + '\t' + 'INFO' + '\t' + 'FORMAT' + '\t' + os.path.basename(BIN).split('.')[0].upper() + '\n'

	return fileformat+filedate+source+contigs+end+motif+refrep+refint+allelerefseq+alleledistance+allelerefnum+allelerefint+altdist+allele1num+allele1int+allele2num+allele2int+gt+gl+dp+ad+head


def VCFV(keyR,REF,ALT,MOTIF,RSIM,RN,RI,RALS,RALD,RALN,RALI,AL1N,AL1I,AL2N,AL2I,GT,GL,DP,AD):

	'''
	Write VCF variant
	'''

	CHROM=keyR.split(':')[0]
	POS=keyR.split(':')[1].split('-')[0]
	ID='.'
	QUAL='.'
	FILTER='.'
	INFO='END=' + keyR.split(':')[1].split('-')[1] + ';MOTIF=' + MOTIF +';RN='+ str(RN) + ';RI='+ RI + ';RALS=' + RALS + ';RALD='+ str(RALD) + ';RALN=' + str(RALN) + ';RALI=' + RALI + ';ALTD='+ str(RSIM) + ';AL1N='+ str(AL1N) + ';AL1I='+ AL1I + ';AL2N='+ str(AL2N) + ';AL2I='+ AL2I
	FORMAT='GT:GL:DP:AD\t' + GT + ':' + GL + ':' + DP + ':' + AD
	variant=CHROM + '\t' + POS + '\t' + ID +'\t' + REF + '\t' + ALT + '\t' + QUAL +'\t' + FILTER +'\t'+ INFO +'\t'+ FORMAT + '\n' 

	return variant


def GTLH(alleles,coverage,error,PHom,PHet):

	'''
	Calculate log10-scaled GT likelihoods for consensus sequences. We do not have per-base quality of the consensus, thus we use the mean erro rate of each group to compute the error probability.
	'''

	listGT=list()
	dictGL=dict()

	for combos in combinations_with_replacement(alleles.keys(),2):

		if combos[0] == combos[1]:

			pC=(PHom**(alleles[combos[0]][2]/coverage))*(error**(sum(y[2] for x,y in alleles.items() if combos[0] != x)/coverage))

		else:

			pC=(PHet**((alleles[combos[0]][2]+alleles[combos[1]][2])/coverage))*(error**(sum(y[2] for x,y in alleles.items() if combos[0] != x and combos[1] != x)/coverage))

		listGT.append(pC)
		dictGL['/'.join(x for x in combos)]=0

	for i,key in enumerate(dictGL.keys()):

		dictGL[key]=np.log10(listGT[i]/sum(x for x in listGT))

	return dictGL



def ParseGroups(BIN,OUT,match,mismatch,gapopen,gapextend,treshold):

	'''
	Generate POA-based consensus sequences for each input group and identify REF/ALT alleles
	'''

	now=datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
	binin=open(BIN,'rb')
	dictR = pickle.load(binin)
	binin.close()
	ctgs=set(x.split(':')[0] for x in list(dictR.keys()))

	#write header and append the other lines afterwards

	with open(os.path.abspath(OUT + '/TREADMILL.vcf'), 'w') as vcfout:

		vcfout.write(VCFH(ctgs,BIN))

	for keyR in dictR.keys():

		now=datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
		print('[' + now + ']' + '[Message] Processing region ' + keyR)

		sMotif=dictR[keyR]['motif']
		OUTR=os.path.abspath(OUT + '/' + keyR)
		dictA=dict()
		listA=list()

		if not os.path.exists(OUTR):

			os.makedirs(OUTR)

		refsequence=dictR[keyR]['reference']
		Rref,Iref=RegExF(refsequence,sMotif)

		with open(os.path.abspath(OUTR + '/r.tmp.fa'), 'w') as fr:

			fr.write('>reference\n' + refsequence + '\n')

		for i,keyG in enumerate(dictR[keyR].keys()):

			if keyG != 'reference' and keyG != 'error' and keyG != 'coverage':

				dictR[keyR][keyG] = OrderedDict(sorted(dictR[keyR][keyG].items(), key=lambda x:x[1][1])) #sort by quality. If the number of sequences is > 100, this allows to exclude low-quality sequences from consensus computation, as these won't be added in the final .fa
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

							listA.append((line.rstrip(), 100*(1-editdistance.eval(line.rstrip(),refsequence)/max(len(line.rstrip()),len(refsequence))), counter))

				os.remove(os.path.abspath(OUTR + '/a'+str(i+1)+'.cs.fa')) #clean-up

		sortlistA=sorted(listA, key=lambda x:x[1], reverse=True) #first is the most similar to the reference

		if sortlistA[0][1] >= treshold:

			RSIM='.'
			RALN,RALI=RegExF(sortlistA[0][0],sMotif)

			for i in range(len(sortlistA)):

				dictA[str(i)] = sortlistA[i]

		else: #if not, all alternative alleles

			RSIM=sortlistA[0][1]
			RALN,RALI='.','.'
			dictA['0'] = ('.','.',0)

			for i in range(len(sortlistA)):

				dictA[str(i+1)] = sortlistA[i]

		#assuming diploidy.

		PHom=1-(len(sortlistA)*dictR[keyR]['error'])
		PHet=.5-dictR[keyR]['error']
		GL=OrderedDict(sorted(GTLH(dictA,dictR[keyR]['coverage'],dictR[keyR]['error'],PHom,PHet).items(),key=lambda x:x[1], reverse=True))
		getKey=next(iter(GL))#get most likely genotype

		if getKey.split('/')[0] != getKey.split('/')[1]:

			if '0' in getKey: #e.g. 0/1, 0/2, 0/3

				getCombos=['/'.join(x) for x in combinations_with_replacement(['0', getKey.split('/')[1]],2)]
				Wgenotype='0/1'
				altal=dictA[getKey.split('/')[1]][0]
				AL1N,AL1I=RegExF(altal,sMotif)
				AL2N,AL2I='.','.'

			else: #e.g. 1/2, 1/3, 2/3, then is multiallelic

				getCombos=list(['0/0', '0/'+getKey.split('/')[0], getKey.split('/')[0]+'/'+getKey.split('/')[0],'0/'+getKey.split('/')[1], getKey, getKey.split('/')[1] + '/' + getKey.split('/')[1]])
				Wgenotype='1/2'
				altal=','.join([dictA[getKey.split('/')[0]][0],dictA[getKey.split('/')[1]][0]])
				AL1N,AL1I=RegExF(dictA[getKey.split('/')[0]][0],sMotif)
				AL2N,AL2I=RegExF(dictA[getKey.split('/')[1]][0],sMotif)

		else: #e.g. 1/1, 2/2

			if getKey.split('/')[0] != '0': #e.g. 1/1, 2/2

				getCombos=['0/0', '0/'+getKey.split('/')[0], getKey.split('/')[0]+'/'+getKey.split('/')[0]]
				Wgenotype='1/1'
				altal=dictA[getKey.split('/')[0]][0]
				AL1N,AL1I=RegExF(altal,sMotif)
				AL2N,AL2I='.','.'

			else: #0/0

				Hkey=[x for x in GL.keys() if '0' in str(x) and x!= '0/0'][0]
				getCombos=['0/0', Hkey,  Hkey.split('/')[1]+'/'+ Hkey.split('/')[1]]
				Wgenotype='0/0'
				altal=''
				AL1N,AL1I='.','.'
				AL2N,AL2I='.','.'

		#store read names in .tsv.gz

		WGL=','.join(str(GL[x]) for x in getCombos)
		WAD=','.join([str(dictA[getKey.split('/')[0]][2]),str(dictA[getKey.split('/')[1]][2])]) 

		originalall1='group' + str(listA.index(dictA[getKey.split('/')[0]])+1) 
		originalall2='group' + str(listA.index(dictA[getKey.split('/')[1]])+1)

		names1=list(dictR[keyR][originalall1].keys())
		names2=list(dictR[keyR][originalall2].keys())

		allnames=names1+names2
		allpos=['1']*len(names1) + ['2']*len(names2)
		allreg=[str(keyR)]*len(allpos)

		now=datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
		print('[' + now + ']' + '[Message] Writing VCF entry')

		with open(os.path.abspath(OUT + '/TREADMILL.vcf'), 'a') as vcfout:

			vcfout.write(VCFV(keyR,refsequence,altal,sMotif,RSIM,Rref,Iref,dictA['0'][0],dictA['0'][1],RALN,RALI,AL1N,AL1I,AL2N,AL2I,Wgenotype,WGL,str(dictR[keyR]['coverage']),WAD))

		now=datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
		print('[' + now + ']' + '[Message] Writing names to TSV')

		tsv=pd.DataFrame({'Region':allreg, 'Name':allnames,'Allele':allpos})
		tsv.to_csv(os.path.abspath(OUT + '/TREADMILL.tsv.gz') ,sep='\t',index=False, header=False, compression='gzip', mode='a')

		shutil.rmtree(OUTR)

	#index

	now=datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
	print('[' + now + ']' + '[Message] Indexing VCF')
	
	pysam.tabix_index(OUT + '/TREADMILL.vcf', preset='vcf', force=True)


def run(parser,args):

	'''
	Execute the code and write variants to VCF.GZ
	'''

	now=datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")

	BIN=os.path.abspath(args.input)

	if not os.path.isfile(BIN):

		print('[' + now + ']' + '[Error] Invalid BIN file')
		sys.exit(1)

	OUT=os.path.abspath(args.output)

	if not os.path.exists(OUT):

		try:

			os.makedirs(OUT)

		except:

			print('[' + now + ']' + '[Error] Cannot create the output folder')
			sys.exit(1)

	ParseGroups(BIN,OUT,args.match,args.mismatch,args.gapopen,args.gapextend,args.similarity)

	now=datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
	print('[' + now + ']' + '[Message] Done')

	sys.exit(0)
