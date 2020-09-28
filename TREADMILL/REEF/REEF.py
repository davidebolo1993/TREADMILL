#!/usr/bin/env python3

import os
import sys
import re
import multiprocessing
import math
import pickle
from itertools import groupby

#additional modules

import pyfaidx
import pysam
import pybedtools
import numpy as np
import mappy as mp


class AutoVivification(dict):

	'''
	Fast way to create nested dictionaries
	'''

	def __getitem__(self, item):
		
		try:
			
			return dict.__getitem__(self, item)
		
		except KeyError:
			
			value = self[item] = type(self)()
			
			return value


def subnone(coordinates):

	'''
	Substitute None (soft-clipped/inserted) coordinates with a negative value
	'''

	return [-9999999 if v is None else v for v in coordinates]


def find_nearest(array, value):

	'''
	Find closest match in array, given value
	'''

	return (np.abs(array - value)).argmin()


def Chunks(l,n):

	'''
	Split list in sublists
	'''

	return [l[i:i+n] for i in range(0, len(l), n)]


def PyCoord(CIGOP,ref_pointer):


	'''
	Extract coords from cigar operations, in the style of pysam (i.e. read.get_reference_positions)
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


def Map(a_instance,map_dict,sequences,flank):

	'''
	Map a list of reads to fake reference sequences
	'''

	for d in sequences:

		seq=d['seq']
		qual=d['qual']
		name=d['name']

		try:

			hit=next(a_instance.map(seq))

			if hit.is_primary:

				w_st=flank #rep start at 500 bp
				w_en=hit.ctg_len-flank #rep end at reference length - flank size

				clip = ['' if x == 0 else '{}S'.format(x) for x in (hit.q_st, len(seq) - hit.q_en)] #calculate soft clipped bases
				cigstr = ''.join((clip[0], hit.cigar_str, clip[1])) #convert to cigarstring	with soft-clipped
				cig_split = [''.join(x) for _, x in groupby(cigstr, key=str.isdigit)] #split by operation
				cigop=[cig_split[n:n+2] for n in range(0,len(cig_split),2)] #group by operation
				cigop_conv=''.join([str(x[1])*int(x[0]) for x in cigop]) #convert to string, one char for each operation
				coord=np.asarray(subnone(PyCoord(cigop_conv,hit.r_st))) #convert to list of coords, in the style of pysam
				si,ei=find_nearest(coord,w_st),find_nearest(coord,w_en)
				subsequence=seq[si:ei]
				suberror=10**(-(np.mean(qual[si:ei]))/10)
				l = map_dict[hit.ctg]
				l.append((name,subsequence,suberror))
				map_dict[hit.ctg]=l
				
			#else:

				#l = map_dict['notprimary']
				#l.append(name)
				#map_dict[hit.ctg]=l

		except StopIteration: #unmapped sequence

			#l = map_dict['unmapped']
			#l.append(name)
			#map_dict[hit.ctg]=l
			continue


def ReMap(BAM,REF,BED,BIN,motifs,flank,maxsize,cores,similarity):

	'''
	Create synthetic chromosomes harboring different set of repeat expansions and map original sequences to these chromosomes
	'''

	hierarchy=AutoVivification()
	fastafile=pyfaidx.Fasta(REF)
	bamfile=pysam.AlignmentFile(BAM, 'rb')
	bedfile=pybedtools.BedTool(BED)
	
	try:
		
		bedsrtd=bedfile.sort()

	except:

		print('[Error] Invalid BED file format')
		sys.exit(1)

	if len(motifs) != len(bedsrtd):

		print('[Error] The number of repeated motifs does not match the number of the regions in the BED file.')
		sys.exit(1)	

	for i,query in enumerate(bedsrtd):

		CHROM,START,END,repeat=query.chrom, query.start,query.end,motifs[i]
		REGION=CHROM + ':' + str(START) + '-' + str(END)

		print('[Message] Processing ' + REGION + ', with repeat ' + motifs[i])

		refseq=fastafile[CHROM][:len(fastafile[CHROM])].seq[START-1:END] #region containing the repetition
		leftflank=fastafile[CHROM][:len(fastafile[CHROM])].seq[START-(flank+1):END-1] #region flanking on the left side
		rightflank=fastafile[CHROM][:len(fastafile[CHROM])].seq[START:END+flank] #region flanking on the right side

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

			REFOUT=os.path.abspath(os.path.dirname(BIN) + '/fake.tmp.fa')

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

			#Adict['unmapped'] = []
			#Adict['notprimary'] = []

			processes=[]
			a=mp.Aligner(REFOUT, preset='map-ont')

			for s in slices:

				p=multiprocessing.Process(target=Map, args=(a,Adict,s,flank))
				p.start()
				processes.append(p)

			for p in processes:
				
					p.join()

			os.remove(REFOUT)

			allerrors=[]
			coverage=0

			for i,key in enumerate(Adict.keys()):

				if not Adict[key] == []: #group not empty

					for el in Adict[key]:

						coverage+=1
						hierarchy[REGION]['group' + str(i+1)][el[0]] = (el[1],el[2])
						allerrors.append(el[2])

			hierarchy[REGION]['reference'] = refseq
			hierarchy[REGION]['coverage'] = coverage
			hierarchy[REGION]['error'] = np.mean(allerrors)

	return hierarchy

		
def run(parser,args):

	'''
	Execute the code and and dump binary output to file
	'''

	BAM=os.path.abspath(args.bamfile)

	if not os.path.isfile(BAM):

		print('[Error] Invalid BAM file')
		sys.exit(1)

	if not os.path.isfile(BAM+'.bai'):

		print('[Warning] Missing BAM file index. Indexing.')
		pysam.index(BAM)

	BED=os.path.abspath(args.bedfile)

	if not os.path.isfile(BED):

		print('[Error] Invalid BED file')
		sys.exit(1)

	REF=os.path.abspath(args.fastafile)

	try:

		with open(REF) as genome:

			assert(genome.readline().startswith('>'))

	except:

		print('[Error] Invalid reference FASTA file')
		sys.exit(1)

	BIN=os.path.abspath(args.output)

	if not os.access(os.path.dirname(BIN),os.W_OK):

		print('[Error] Missing write permissions on the output folder')
		sys.exit(1)

	hierarchy=ReMap(BAM,REF,BED,BIN,args.motif[0],args.flanking,args.maxsize,args.threads,args.similarity)

	binout=open(BIN,'wb')
	data=pickle.dumps(hierarchy,protocol=pickle.HIGHEST_PROTOCOL)
	binout.write(data)
	binout.close()

	#reopening for checking routine
	#binin=open(args.output,'rb')
	#data = pickle.load(binin)
	#binin.close()

	sys.exit(0)

