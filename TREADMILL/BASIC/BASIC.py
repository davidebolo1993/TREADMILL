#!/usr/bin/env python3

#python 3 standard library

import sys
import os
import json
import gzip
from collections import defaultdict
from datetime import datetime

#additional modules

import pysam
import pybedtools
import numpy as np

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

	S_dict=dict()

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
	S_dict['BAM_LEN'] = [] #all lengths in list
	S_dict['BAM_QUAL'] = [] #all qualities in list
	S_dict['BAM_PID'] = [] #all PID in list

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')

	bamfile=pysam.AlignmentFile(BAM, 'rb')
	bedfile=pybedtools.BedTool(BED)
	
	try:
		
		bedsrtd=bedfile.sort()

	except:

		print('[' + now + ']' + '[Error] Invalid BED file format')
		sys.exit(1)

	print('[' + now + ']' + '[Message] Parsing BAM file')

	ivf=bedsrtd.as_intervalfile()

	for read in bamfile.fetch():

		if read.has_tag('MD') and read.has_tag('NM'):

			if not read.is_unmapped:

				if not read.is_supplementary:

					if not read.is_secondary:

						S_dict['BAM_PRIM']+=1

						query=pybedtools.Interval(read.reference_name,read.reference_start,read.reference_end)
						
						if ivf.any_hits(query) >=1: #parse CIGAR only in targeted regions

							S_dict['BAM_ONTARGET']+=1

							Cdict=tuptodict(read.cigartuples)
							MD=read.get_aligned_pairs(with_seq=True)

							S_dict['BAM_CSOFT_CLIP']+=Cdict[4]
							S_dict['BAM_CINS']+=Cdict[1]
							S_dict['BAM_CDEL']+=Cdict[2]
							S_dict['BAM_CMATCH']+=sum(1 for x,y,z in MD if x is not None and z is not None and z[0].isupper())
							S_dict['BAM_CDIFF']+=sum(1 for x,y,z in MD if x is not None and z is not None and z[0].islower())
							S_dict['BAM_QUAL'].append(np.mean(read.query_qualities))

							NM=read.get_tag('NM')
							refcoords=read.get_reference_positions()
							reflen=refcoords[-1]-refcoords[0]#reference spans from first aligned to last aligned
							seqlen=len(read.get_reference_positions(full_length=True))#this is the read length
							S_dict['BAM_LEN'].append(seqlen)
							PID=100-100*NM/max(reflen,seqlen)
							S_dict['BAM_PID'].append(PID)

						else:

							S_dict['BAM_OFFTARGET']+=1

					else:

						S_dict['BAM_SEC']+=1

				else:

					S_dict['BAM_SUPP']+=1

			else:

				S_dict['BAM_UNMAP']+=1

		else:

			print('[Error] BAM misses the required MD/NM tags')
			bamfile.close()
			sys.exit(1)

	#calculate coverage in regions from BED

	for query in bedsrtd:

		key=query.chrom+':'+str(query.start)+'-'+str(query.end)
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Message] Calculating coverage in region ' + key)
		query_arr=bamfile.count_coverage(query.chrom,query.start,query.end,read_callback=check_read)
		perbasecov=np.sum(query_arr,axis=0).tolist()
		S_dict[key]=perbasecov

	bamfile.close()

	return S_dict


def run(parser,args):

	'''
	Execute the code and dump JSON to file
	'''

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] TREADMILL BASIC v1.0')

	BAM=os.path.abspath(args.bamfile)

	if not os.path.isfile(BAM):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid BAM file')
		sys.exit(1)

	if not os.path.isfile(BAM+'.bai'):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Missing BAM file index')
		sys.exit(1)

	BED=os.path.abspath(args.bedfile)

	if not os.path.isfile(BED):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid BED file')
		sys.exit(1)

	JSON=os.path.abspath(args.output)

	if not os.access(os.path.dirname(JSON),os.W_OK):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Missing write permissions on the output folder')
		sys.exit(1)


	S_dict=parse(BAM,BED)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Writing output')

	with gzip.open(JSON+'.gz', 'wt') as gzout:

		json.dump(S_dict, gzout, indent=4)
		gzout.write('\n')

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Done')

	sys.exit(0)