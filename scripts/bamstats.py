#!/usr/bin/env python3

#python 3 standard library

import sys
import os
import re
import json
from collections import defaultdict

#additional modules

import pysam
import pybedtools
from pybedtools import Interval
import numpy as np

#global

S_dict=dict()

def tuptodict(cigartuples):

	'''
	Convert tuples from read.cigartuples to dict and sum values having same keys
	'''

	Cdict = defaultdict(int)

	for key, value in cigartuples:

		Cdict[key] = Cdict[key]+value

	return Cdict



def check_read(read):

	'''
	Return True for primary alignments
	'''

	if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

		return True


def parse(BAM,BED): 

	'''
	Parse BAM and get overall and queries-specific statistics
	'''

	global S_dict

	S_dict['BAM_PRIM']=0 #reads
	S_dict['BAM_SUPP']=0 #reads
	S_dict['BAM_SEC']=0 #reads
	S_dict['BAM_UNMAP']=0 #reads
	S_dict['BAM_ONTARGET']=0 #reads
	S_dict['BAM_OFFTARGET']=0 #reads
	S_dict['BAM_CSOFT_CLIP']=0 #bps
	S_dict['BAM_CMATCH']=0 #bps
	S_dict['BAM_CINS']=0 #bps
	S_dict['BAM_CDEL']=0 #bps
	S_dict['BAM_CDIFF']=0 #bps

	bamfile=pysam.AlignmentFile(BAM, "rb")
	bed=pybedtools.BedTool(BED)
	bedsrtd=bed.sort()
	ivf=bedsrtd.as_intervalfile()

	for read in bamfile.fetch():

		if read.has_tag('MD'):

			if not read.is_unmapped:

				if not read.is_supplementary:

					if not read.is_secondary:

						S_dict['BAM_PRIM']+=1

						query=Interval(read.reference_name,read.reference_start,read.reference_end)
						
						if ivf.any_hits(query) >=1: #parse CIGAR only in targeted regions

							S_dict['BAM_ONTARGET']+=1

							Cdict=tuptodict(read.cigartuples)
							MD=read.get_aligned_pairs(with_seq=True)

							S_dict['BAM_CSOFT_CLIP']+=Cdict[4]
							S_dict['BAM_CINS']+=Cdict[1]
							S_dict['BAM_CDEL']+=Cdict[2]
							S_dict['BAM_CMATCH']+=sum(1 for x,y,z in MD if x is not None and z is not None and z[0].isupper())
							S_dict['BAM_CDIFF']+=sum(1 for x,y,z in MD if x is not None and z is not None and z[0].islower())

						else:

							S_dict['BAM_OFFTARGET']+=1

					else:

						S_dict['BAM_SEC']+=1

				else:

					S_dict['BAM_SUPP']+=1

			else:

				S_dict['BAM_UNMAP']+=1

		else:

			print('[Error] Missing required MD tag')
			bamfile.close()
			sys.exit(1)


	#calculate coverage in regions from BED

	for query in bedsrtd:

		query_arr=bamfile.count_coverage(query.chrom,query.start,query.end,read_callback=check_read)
		meancov=np.mean(np.sum(query_arr,axis=0))
		key=query.chrom+":"+str(query.start)+"-"+str(query.end)
		S_dict[key]=meancov


	bamfile.close()

	return S_dict


def main():

	'''
	Execute the code and dump json to stdout
	'''

	try:
		
		assert(len(sys.argv) == 3)

	except:

		print("[Usage] python bamstats.py <input.bam> <queries.bed>")
		sys.exit(1)

	BAM=os.path.abspath(sys.argv[1])

	if not os.path.isfile(BAM):

		print('[Error] Invalid BAM file <input.bam>')
		print("[Usage] python bamstats.py <input.bam> <queries.bed>")
		sys.exit(1)


	BED=os.path.abspath(sys.argv[2])

	if not os.path.isfile(BED):

		print('[Error] Invalid BED file <queries.bed>')
		print("[Usage] python bamstats.py <input.bam> <queries.bed>")
		sys.exit(1)

	S_dict=parse(BAM,BED)
	json.dump(S_dict, sys.stdout)

	sys.exit(0)

if __name__=='__main__':

	main()
