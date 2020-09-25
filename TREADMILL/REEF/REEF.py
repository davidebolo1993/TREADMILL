#!/usr/bin/env python3

import os
import sys
import re
import multiprocessing
import math
from itertools import groupby

#additional modules

import pyfaidx
import pysam
import numpy as np
import mappy as mp




def subnone(coordinates):

	'''
	Substitute None (soft-clipped/inserted) coordinates with a negative value
	'''

	return [-999999 if v is None else v for v in coordinates]


def find_nearest(array, value):

	'''
	Find closest match in array, given value
	'''

	return (np.abs(array - value)).argmin()



def Chunks(l,n):

	'''
	Split list l in n sublists
	'''

	return [l[i:i+n] for i in range(0, len(l), n)]




def PyCoord(CIGOP,ref_pointer):


	'''
	Extract coords from cigar operations, in the style of pysam
	'''

	s=ref_pointer

	coords=[]

	for op in list(CIGOP):

		if op == 'M':

			coords.append(s)
			s+=1

		elif op == 'I' or op == 'S':

			coords.append(None)

		elif op == 'D':

			s+=1

	return coords

		



def Map(a_instance,map_dict,sequences):

	'''
	Map a list of sequences to a reference index
	'''

	for d in sequences:

		seq=d['seq']
		qual=d['qual']
		name=d['name']

		try:

			hit=next(a_instance.map(seq))

			if hit.is_primary:

				w_st=500 #rep start at 500 bp
				w_en=hit.ctg_len-500 #rep end at reference length - 500 bp

				clip = ['' if x == 0 else '{}S'.format(x) for x in (hit.q_st, len(seq) - hit.q_en)] #calculate soft clipped bases
				cigstr = ''.join((clip[0], hit.cigar_str, clip[1])) #convert to cigarstring	
				cig_split = [''.join(x) for _, x in groupby(cigstr, key=str.isdigit)]
				cigop=[cig_split[n:n+2] for n in range(0,len(cig_split),2)]
				cigop_conv=''.join([str(x[1])*int(x[0]) for x in cigop])
				coord=np.asarray(subnone(PyCoord(cigop_conv,hit.r_st)))
				si,ei=find_nearest(coord,w_st),find_nearest(coord,w_en)
				subsequence=seq[si:ei]
				suberror=10**(-(np.mean(qual[si:ei]))/10)
				l = map_dict[hit.ctg]
				l.append((name,subsequence,suberror))
				map_dict[hit.ctg]=l
			
				print(hit)

			else:

				l = map_dict['notprimary']
				l.append(name)
				map_dict[hit.ctg]=l

		except StopIteration: #unmapped sequence

			l = map_dict['unmapped']
			l.append(name)
			map_dict[hit.ctg]=l


def ReMap(BAM,REF,REGION,OUT,repeat,maxsize,cores,similarity):

	'''

	'''

	fastafile=pyfaidx.Fasta(REF)
	bamfile=pysam.AlignmentFile(BAM, 'rb')

	CHROM,START,END=REGION.split(':')[0],int(REGION.split(':')[1].split('-')[0]),int(REGION.split(':')[1].split('-')[1])

	refseq=fastafile[CHROM][:len(fastafile[CHROM])].seq[START-1:END] #region containing the repetition strictly
	leftflank=fastafile[CHROM][:len(fastafile[CHROM])].seq[START-501:END-1] #500 bps region flanking on the left side
	rightflank=fastafile[CHROM][:len(fastafile[CHROM])].seq[START:END+500] #500 bps region flanking on the right side

	r=re.compile(r'('  + repeat + r')\1+')
	seen=set()

	for match in r.finditer(refseq):

		seen.add(match.group(0))

	if len(seen) == 0:

		print('[Error] Could not find repeat in reference')
		sys.exit(1)

	else:

		min_=sum(el.count(repeat) for el in seen)
		print('[Message] ' + str(min_) + ' ' + str(repeat) + ' in reference sequence')

		seqdict=dict()

		header='>treadmill_reef_' + REGION + '_' + repeat + '_' + str(min_)
		sequence=leftflank+refseq+rightflank
		seqdict[header]=sequence
		reps=min_

		while maxsize - reps > min_:

			seqtoadd=(len(sequence)/100)*(100-similarity)
			reptoadd=round(seqtoadd/len(repeat))
			reps+=reptoadd
			header='>treadmill_reef_' + REGION + '_' + repeat + '_' + str(reps)
			sequence=leftflank+refseq+(repeat*(reps-min_))+rightflank
			seqdict[header]=sequence

		REFOUT=os.path.abspath(OUT + '/fake.fa')

		with open(REFOUT, 'w') as fout:

			for key,value in seqdict.items():

				fout.write(key+'\n'+value+'\n')

		BAMseqs=list()

		for reads in bamfile.fetch(CHROM,START,END):

			if not reads.is_unmapped and not reads.is_supplementary and not reads.is_secondary: #just primary aligments

				Rdict=dict()

				Rdict['name'] = reads.query_name
				Rdict['seq'] = reads.query_sequence
				Rdict['qual'] = reads.query_qualities

				BAMseqs.append(Rdict)

		chunk_size=len(BAMseqs)/cores
		slices=Chunks(BAMseqs,math.ceil(chunk_size))
		manager=multiprocessing.Manager()
		Adict=manager.dict()

		for key in seqdict.keys(): #intialize empty

			Adict[key[1:]] = []

		Adict['unmapped'] = []
		Adict['notprimary'] = []

		processes=[]
		a=mp.Aligner(REFOUT, preset='map-ont')

		for s in slices:

			p=multiprocessing.Process(target=Map, args=(a,Adict,s))
			p.start()
			processes.append(p)

		for p in processes:
			
				p.join()

		




def run(parser,args):

	'''
	Execute the code and ....
	'''

	BAM=os.path.abspath(args.bamfile)

	if not os.path.isfile(BAM):

		print('[Error] Invalid BAM file')
		sys.exit(1)

	if not os.path.isfile(BAM+'.bai'):

		print('[Warning] Missing BAM file index. Indexing.')
		pysam.index(BAM)


	REF=os.path.abspath(args.fastafile)

	try:

		with open(REF) as genome:

			assert(genome.readline().startswith('>'))

	except:

		print('[Error] Invalid reference FASTA file')
		sys.exit(1)

	r=re.compile('.*:.*-.*')

	if r.match(args.region[0]) is None:

		print('[Error] Wrong region format')
		sys.exit(1)

	OUT=os.path.abspath(args.output)

	if not os.path.exists(OUT):

		os.makedirs(OUT)


	ReMap(BAM,REF,args.region[0],OUT, args.repeat, args.maxsize,args.threads,args.similarity)

