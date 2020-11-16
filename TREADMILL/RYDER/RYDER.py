#!/usr/bin/env python3

import os
import sys
import re
import multiprocessing
import math
import pickle
from itertools import groupby,combinations
from datetime import datetime

#additional modules

import pyfaidx
import pysam
import pybedtools
import numpy as np
import mappy as mp
import editdistance

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


def similarity(worda,wordb):

	'''
	Return the edit distance-based similarity score between 2 sequences
	'''

	return 100-100*editdistance.eval(worda,wordb)/max(len(worda),len(wordb))


def decisiontree(readsdict,mingroupsize,treshold):

	'''
	Cluster strings in list by similarity (edit distance score)
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

	for i,op in enumerate(list(CIGOP)):

		if i==0 and op == 'M':

			coords.append(s)

		else:

			if op == 'M':

				s+=1
				coords.append(s)

			elif op == 'I' or op == 'S':

				coords.append(None)

			elif op == 'D':

				s+=1

	return coords


def BamW(header,BAMsegments,BAM):

	'''
	Write aligned segments in BAM format

	'''

	with pysam.AlignmentFile(BAM, mode='wb', header=header) as bout:

		for segments in BAMsegments:

			s = pysam.AlignedSegment(bout.header)
			s.is_unmapped=False #'cause unmapped reads were skipped
			s.is_reverse=False #'cause reads are all translated to forward orientation
			s.is_secondary=False #'cause only primary alignments were retained
			s.query_name=segments['QNAME']
			s.reference_name=segments['RNAME']
			s.reference_start=segments['POS']
			s.mapping_quality=segments['MAPQ']
			s.cigarstring=segments['CIGAR']
			s.query_sequence=segments['SEQ']
			s.query_qualities=segments['QUAL']
			s.set_tags([('MD',segments['MD'], 'Z'), ('cs', segments['cs'], 'Z')])
			bout.write(s)


def Map(a_instance,Slist,Qlist,sequences,flank,finalBAM,store):

	'''
	Map a list of reads to fake reference sequences
	'''

	for d in sequences:

		seq=d['seq']
		qual=d['qual']
		name=d['name']

		Aldict=dict()

		try:

			hit=next(a_instance.map(seq,MD=True,cs=True))

			if hit.is_primary:

				w_st=flank #rep start at flank
				w_en=hit.ctg_len-flank #rep end at reference length - flank size

				clip = ['' if x == 0 else '{}S'.format(x) for x in (hit.q_st, len(seq) - hit.q_en)] #calculate soft clipped bases
				cigstr = ''.join((clip[0], hit.cigar_str, clip[1])) #convert to cigarstring	with soft-clipped
				cig_split = [''.join(x) for _, x in groupby(cigstr, key=str.isdigit)] #split by operation
				cigop=[cig_split[n:n+2] for n in range(0,len(cig_split),2)] #group by operation
				cigop_conv=''.join([str(x[1])*int(x[0]) for x in cigop]) #convert to string, one char for each operation
				coord=np.asarray(subnone(PyCoord(cigop_conv,hit.r_st))) #convert to list of coords, in the style of pysam
				si,ei=find_nearest(coord,w_st),find_nearest(coord,w_en) 

				if hit.r_st <= w_st and hit.r_en >= w_en: #the entire rep is spanned, keep the read

					subsequence=seq[si:ei]
					suberror=10**(-(np.mean(qual[si:ei]))/10)
					
					Slist.append((name,subsequence))
					Qlist.append((name, suberror))

				if store:

					Aldict['QNAME'] = name
					Aldict['RNAME'] = hit.ctg
					Aldict['POS'] = hit.r_st				
					Aldict['MAPQ'] = hit.mapq
					Aldict['CIGAR'] = cigstr
					Aldict['SEQ'] = seq
					Aldict['QUAL'] = qual
					Aldict['MD'] = hit.MD
					Aldict['cs'] = hit.cs

					finalBAM.append(Aldict)
				
			else:

				pass
				
				
		except StopIteration: #unmapped sequence
			
			pass


def ReMap(BAM,REF,BED,BIN,motifs,flank,maxsize,cores,sim,support,store):

	'''
	Create synthetic chromosomes harboring different set of repeat expansions and map original sequences to these chromosomes. Group reads by similarity.
	'''

	hierarchy=AutoVivification()
	fastafile=pyfaidx.Fasta(REF)
	bamfile=pysam.AlignmentFile(BAM, 'rb')
	
	try:
		
		bedfile=pybedtools.BedTool(BED)
		bedsrtd=bedfile.sort()

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid BED file format')
		sys.exit(1)

	FAKEREF=os.path.abspath(os.path.dirname(BIN) + '/fake.fa')
	FAKEBAM=os.path.abspath(os.path.dirname(BIN) + '/fake.bam')
	FAKESRTBAM=os.path.abspath(os.path.dirname(BIN) + '/fake.srt.bam')

	manager=multiprocessing.Manager()
	BAMsegments=manager.list() #initialize even if this will stay empty
	
	if len(motifs) != len(bedsrtd):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] The number of repeated motifs provided does not match the number of regions in the BED file')
		sys.exit(1)	

	for i,query in enumerate(bedsrtd):

		CHROM,START,END,repeat=query.chrom, query.start,query.end,motifs[i]
		REGION=CHROM + ':' + str(START) + '-' + str(END)
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Message] Processing region ' + REGION+ ', with repeat ' + motifs[i])

		refseq=fastafile[CHROM][:len(fastafile[CHROM])].seq[START-1:END] #region containing the repetition
		leftflank=fastafile[CHROM][:len(fastafile[CHROM])].seq[START-(flank+1):END-1] #region flanking on the left side
		rightflank=fastafile[CHROM][:len(fastafile[CHROM])].seq[START:END+flank] #region flanking on the right side

		seen=[m.start() for m in re.finditer(repeat, refseq)]
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')

		if len(seen) == 0:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Error] Could not find the specified repeat in reference')
			sys.exit(1)

		else:

			min_=len(seen)
			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Message] ' + str(min_) + ' ' + repeat + ' in reference sequence')

			seqdict=dict()

			header='>treadmill_reef_' + REGION + '_' + repeat + '_' + str(min_) + '_flank_' + str(flank)
			sequence=leftflank+refseq+rightflank
			seqdict[header]=sequence
			reps=min_

			while maxsize - reps > min_:

				seqtoadd=(len(sequence)/100)*(100-sim)
				reptoadd=round(seqtoadd/len(repeat))
				reps+=reptoadd
				header='>treadmill_reef_' + REGION + '_' + repeat + '_' + str(reps) + '_flank_' + str(flank)
				sequence=leftflank+refseq+(repeat*(reps-min_))+rightflank
				seqdict[header]=sequence

			REFOUT=os.path.abspath(os.path.dirname(BIN) + '/fake.tmp.fa')

			with open(REFOUT, 'w') as fout:

				for key,value in seqdict.items():

					fout.write(key+'\n'+value+'\n')

			BAMseqs=list()
			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Message] Parsing BAM file')

			for reads in bamfile.fetch(CHROM,START,END):

				if not reads.is_unmapped and not reads.is_supplementary and not reads.is_secondary: #just primary aligments

					Rdict=dict()

					Rdict['name'] = reads.query_name
					Rdict['seq'] = reads.query_sequence
					Rdict['qual'] = reads.query_qualities

					BAMseqs.append(Rdict)

			#parallelize alignment

			chunk_size=len(BAMseqs)/cores
			slices=Chunks(BAMseqs,math.ceil(chunk_size))
			Slist=manager.list()
			Qlist=manager.list()

			processes=[]
			a=mp.Aligner(REFOUT, preset='map-ont')

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Message] Re-mapping reads to the synthetic reference')

			for s in slices:

				p=multiprocessing.Process(target=Map, args=(a,Slist,Qlist,s,flank,BAMsegments,store))
				p.start()
				processes.append(p)

			for p in processes:
				
				p.join()

			#append to fake reference if store is True and clean-up

			if store:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + ']' + '[Message] Writing synthetic contigs to file')

				faref = open(FAKEREF, 'at')
				txref = open(REFOUT, 'rt')
				faref.write(txref.read()) 			
				faref.close()
				txref.close()

			os.remove(REFOUT)

			#convert list of tuples to dict for better usability

			sdict=dict((x, y) for x, y in Slist)
			qdict=dict((x, y) for x, y in Qlist)

			#fine tuning: re-group by similarity of sequences. This avoids having outlayers that decrease consensus accuracy in TRAP.

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Message] Grouping reads by similarity')

			decision=decisiontree(sdict,support,sim)
			allerrors=[]

			for i,groups in enumerate(decision):

				group='group'+str(i+1)

				for elements in groups:

					keys=[k for k,v in sdict.items() if v == elements] #get corresponding keys
				
					for k in keys:

						hierarchy[REGION][group][k]=(sdict[k],qdict[k])
						allerrors.append(qdict[k])

			hierarchy[REGION]['reference'] = refseq
			hierarchy[REGION]['coverage'] = len(sdict)
			hierarchy[REGION]['error'] = np.mean(allerrors)
			hierarchy[REGION]['motif'] = repeat

	if store: #also write to BAM if store is True

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Message] Writing re-aligned segments to file')
		faref=pyfaidx.Fasta(FAKEREF)
		header = {'HD': {'VN': 1.6, 'SO': 'coordinate'},'SQ': [{'LN': len(faref[k]),'SN': k} for k in faref.keys()]}
		BamW(header,BAMsegments,FAKEBAM)
		pysam.sort('-o', FAKESRTBAM, '-@', str(cores), FAKEBAM)
		pysam.index(FAKESRTBAM)
		os.remove(FAKEBAM)

	return hierarchy

		
def run(parser,args):

	'''
	Execute the code and and dump binary output to file
	'''

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] TREADMILL RYDER v1.0')

	if not os.path.exists(BAM):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid BAM file')
		sys.exit(1)

	if not os.path.exists(BAM+'.bai'):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Warning] Missing BAM file index. Creating')
		
		try:

			pysam.index(BAM)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Error] Cannot index BAM')
			sys.exit(1)

	BED=os.path.abspath(args.bedfile)

	if not os.path.exists(BED):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid BED file')
		sys.exit(1)

	REF=os.path.abspath(args.fastafile)

	try:

		with open(REF) as genome:

			assert(genome.readline().startswith('>'))

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid reference FASTA file')
		sys.exit(1)

	BIN=os.path.abspath(args.output)

	if not os.access(os.path.dirname(BIN),os.W_OK):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Missing write permissions on the output folder')
		sys.exit(1)

	hierarchy=ReMap(BAM,REF,BED,BIN,args.motif[0],args.flanking,args.maxsize,args.threads,args.similarity,args.support,args.store)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Writing output')

	binout=open(BIN,'wb')
	data=pickle.dumps(hierarchy,protocol=pickle.HIGHEST_PROTOCOL)
	binout.write(data)
	binout.close()

	#reopening for checking routine
	#binin=open(args.output,'rb')
	#data = pickle.load(binin)
	#binin.close()

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Done')

	sys.exit(0)