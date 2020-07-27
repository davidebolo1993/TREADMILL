#!/usr/bin/env python3

#python 3 standard library

import sys
import os
import pickle
from itertools import combinations

#additional modules

import pysam
import editdistance
import pybedtools
import numpy as np


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
	Substitute None (soft-clipped/insertions) coordinates with a negative value
	'''

	return [-999999 if v is None else v for v in coordinates]


def find_nearest(array, value):

	'''
	Find closest match in array, given value
	'''

	return (np.abs(array - value)).argmin()



def similarity(worda,wordb):

	'''
	Return the edit-distance based similarity score between 2 sequences
	'''

	return 100-100*editdistance.eval(worda,wordb)/(len(worda)+len(wordb))


def decisiontree(readsdict,mingroupsize,treshold):

	'''
	Group strings in list by similarity (edit distance score)
	'''

	paired ={c:{c} for c in readsdict.values()}

	for worda,wordb in combinations(readsdict.values(),2):

		if similarity(worda,wordb) < treshold: 

				continue

		else:

			paired[worda].add(wordb)
			paired[wordb].add(worda)


	decision = list()
	ungrouped = set(readsdict.values())
	
	while ungrouped:

		best = {}

		for word in ungrouped:

			g = paired[word] & ungrouped

			for c in g.copy():
			
				g &= paired[c]

			if len(g) > len(best):

				best=g

		if len(best) < mingroupsize:

			break

		ungrouped -= best

		decision.append(best)

	return decision



def parser(BAM,BED,mingroupsize,treshold):

	'''
	Parse BAM file and group read by similarity
	'''

	hierarchy=AutoVivification()
	bedfile=pybedtools.BedTool(BED)
	
	try:

		bedsrtd=bedfile.sort() #this allows to also check for errors in the BED

	except:

		print('[Error] Invalid BED file format')
		sys.exit(1)

	bamfile=pysam.AlignmentFile(BAM, 'rb')

	for query in bedsrtd:

		key=query.chrom+':'+str(query.start)+'-'+str(query.end)
		sdict=dict()
		qdict=dict()

		for read in bamfile.fetch(query.chrom,query.start,query.end):

			if not read.is_unmapped and not read.is_supplementary and not read.is_secondary:

				if read.reference_start <= query.start and read.reference_end >= query.end: #retain reads that span entire repetition

					identifier=read.query_name
					sequence=read.query_sequence #no need to reverse complement if read.is_reverse, as query_sequence always return forward
					quality=read.query_qualities
					coord=np.asarray(subnone(read.get_reference_positions(full_length=True)))
					si,ei=find_nearest(coord,query.start),find_nearest(coord,query.end)
					subsequence=sequence[si:ei+1]
					suberror=10**(-(np.mean(quality[si:ei+1]))/10) #extract error probability
					sdict[identifier]=subsequence
					qdict[identifier]=suberror

		decision=decisiontree(sdict,mingroupsize,treshold)

		for i,groups in enumerate(decision):

			group='group'+str(i+1)

			for elements in groups:

				keys=[k for k,v in sdict.items() if v == elements] #get corresponding keys
				
				for k in keys:

					hierarchy[key][group][k][sdict[k]]=qdict[k]

	bamfile.close()

	return hierarchy


def run(parser,args):

	'''
	Execute the code and dump binary output to file
	'''

	BAM=os.path.abspath(args.bamfile)

	if not os.path.isfile(BAM):

		print('[Error] Invalid BAM file')
		sys.exit(1)

	if not os.path.isfile(BAM+'.bai'):

		print('[Error] Missing BAM file index')
		sys.exit(1)

	BED=os.path.abspath(args.bedfile)

	if not os.path.isfile(BED):

		print('[Error] Invalid BED file')
		sys.exit(1)

	BIN=os.path.abspath(args.output)

	if not os.access(os.path.dirname(BIN),os.W_OK):

		print('[Error] Missing write permissions on the output folder')
		sys.exit(1)

	hierarchy=parser(BAM,BED,args.support,args.similarity)

	binout=open(args.output,'wb')
	data=pickle.dumps(hierarchy,protocol=pickle.HIGHEST_PROTOCOL)
	binout.write(data)
	binout.close()

	#reopening for checking routine
	#binin=open(args.output,'rb')
	#data = pickle.load(binin)
	#binin.close()

	sys.exit(0)