#!/usr/bin/env python3

#python 3 standard library

import sys
import os
import json
import gzip
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

			print('[Error] BAM is missing required MD tag')
			bamfile.close()
			sys.exit(1)


	#calculate coverage in regions from BED

	for query in bedsrtd:

		query_arr=bamfile.count_coverage(query.chrom,query.start,query.end,read_callback=check_read)
		perbasecov=np.sum(query_arr,axis=0).tolist()
		key=query.chrom+":"+str(query.start)+"-"+str(query.end)
		S_dict[key]=perbasecov


	bamfile.close()

	return S_dict


def run(parser,args):

	'''
	Execute the code and dump json to stdout
	'''

	BAM=os.path.abspath(args.bamfile)

	if not os.path.isfile(BAM):

		print('[Error] Invalid BAM file')
		sys.exit(1)

	BED=os.path.abspath(args.bedfile)

	if not os.path.isfile(BED):

		print('[Error] Invalid BED file')
		sys.exit(1)

	JSON=os.path.abspath(args.output)

	if not os.access(os.path.dirname(JSON),os.W_OK):

		print('[Error] Missing write permissions on the output folder')
		sys.exit(1)

	S_dict=parse(BAM,BED)

	if args.gzipped:

		with gzip.GzipFile(JSON+".gz", 'w') as gzout:
   
			gzout.write(json.dumps(S_dict).encode('utf-8'))  

	else:

		with open(JSON, 'w') as plainout:

			palinout.write(json.dumps(S_dict).encode('utf-8'))  

	sys.exit(0)