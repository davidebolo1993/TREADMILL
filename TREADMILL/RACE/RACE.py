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
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import AgglomerativeClustering,DBSCAN


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

	return [-9999999999 if v is None else v for v in coordinates]


def find_nearest(array, value):

	'''
	Find closest match in array, given value
	'''

	return (np.abs(array - value)).argmin()


def diffperc(x,y):

	'''
	Custom function for string clustering
	'''

	return 100*int(editdistance.eval(data[int(x[0])], data[int(y[0])]))/max(len(data[int(x[0])]),len(data[int(y[0])]))



def editdist(x,y):

	'''
	Custom function for string clustering
	'''

	return int(editdistance.eval(data[int(x[0])], data[int(y[0])]))


def decisiontree(readsdict,mingroupsize,cluster,tresh):

	'''
	Cluster strings
	'''

	global data
	data=list(readsdict.values())
	X = np.arange(len(data)).reshape(-1, 1)
	result=[]

	if not cluster: #perform clustering based on string similarity

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Message] Performing DBSCAN')
		print('[' + now + ']' + '[Message] Clustering reads by similarity')

		metric= pairwise_distances(X, X, metric=diffperc)
		agg=DBSCAN(eps=100-tresh, min_samples=mingroupsize,algorithm='brute', metric='precomputed')

	else:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Message] Performing agglomerative clustering')

		metric= pairwise_distances(X, X, metric=editdist)

		if type(tresh) == int: #by number of clusters

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Message] Creating ' + str(tresh) + ' clusters')
			agg = AgglomerativeClustering(n_clusters=tresh, affinity='precomputed',linkage='average')

		elif type(tresh) == float: #by threshold

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Message] Cutting dendogram tree using threshold ' + str(tresh))
			agg = AgglomerativeClustering(n_clusters=None,distance_threshold=tresh,affinity='precomputed',linkage='average')

		else: #is string. Compute full dendogram ans perform Silhouette analysis

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')			
			print('[' + now + ']' + '[Message] Computing full dendogram and storing model to output file')
			agg=AgglomerativeClustering(distance_threshold=0, n_clusters=None, affinity='precomputed', linkage='average')

	cluster_=agg.fit(metric)

	if type(tresh) == str: #full tree was computed and model can be plotted

		return [cluster_,metric]

	else:

		groups=set(cluster_.labels_)

		for g in groups:

			if g == -1:

				continue

			else:

				group=list(np.take(data,np.where(cluster_.labels_ == g))[0])

				if len(group) >= mingroupsize: #this only applies to clustering not DBSCAN in practice

					result.append(group)

	return result,metric


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


def ReMap(BAM,REF,BED,BIN,motifs,flank,maxsize,cores,sim,support,store,cluster,tresh):

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

	if cluster:

		if len(bedsrtd) > 1:

			if type(tresh) == str:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + ']' + '[Error] Only one region at a time must be provided in BED when computing the full dendogram')
				sys.exit(1)

			else: #this is int or float

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + ']' + '[Warning] Value given to --clusters/--threshold will be propagated to all the regions in BED')

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
		leftflank=fastafile[CHROM][:len(fastafile[CHROM])].seq[START-(flank+1):START-1] #region flanking on the left side
		rightflank=fastafile[CHROM][:len(fastafile[CHROM])].seq[END:END+flank] #region flanking on the right side

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
			print('[' + now + ']' + '[Message] Clustering reads')

			decision,silhouette=decisiontree(sdict,support,cluster,tresh)

			if type(tresh) == str: #then decision is a model

				hierarchy=decision

			else: #decision is actual the clustered result

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

	return hierarchy,silhouette

		
def run(parser,args):

	'''
	Execute the code and and dump binary output to file
	'''

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] TREADMILL RACE v1.0')

	BAM=os.path.abspath(args.bamfile)

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


	OUTDIR=os.path.dirname(os.path.abspath(args.output))

	if not os.path.exists(OUTDIR):

		try:

			os.makedirs(OUTDIR)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Error] Cannot create the output folder')
			sys.exit(1)

	BIN=os.path.abspath(args.output)

	if not os.access(os.path.dirname(BIN),os.W_OK):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Missing write permissions on the output folder')
		sys.exit(1)

	cluster=False

	if args.hierarchical_clustering:

		if not args.threshold and not args.clusters and not args.dendogram:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Error] When performing hierarchical clustering, one between --threshold, --clusters and --dendogram must be specified')
			sys.exit(1)

		elif (args.threshold and args.clusters) or (args.threshold and args.dendogram) or (args.clusters and args.dendogram):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Error] When performing hierarchical clustering, only one between --threshold, --clusters and --dendogram must be specified')
			sys.exit(1)

		else:

			cluster=True
			
			if args.threshold:

				tresh=float(args.threshold)

			elif args.clusters:

				tresh=int(args.clusters)

			else: #args.dendogram

				tresh='dendogram'
				silout=os.path.abspath(OUTDIR + '/simmatrix.bin')

	else: #greedy string clustering

		tresh=args.affinity

	hierarchy,silhouette=ReMap(BAM,REF,BED,BIN,args.motif[0],args.flanking,args.maxsize,args.threads,args.similarity,args.support,args.store,cluster,tresh)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Writing output')

	binout=open(BIN,'wb')
	data=pickle.dumps(hierarchy,protocol=pickle.HIGHEST_PROTOCOL)
	binout.write(data)
	binout.close()

	if tresh=='dendogram':

		with open(silout, 'wb') as sout:

			np.save(sout,silhouette)

	#reopening for checking routine
	#binin=open(args.output,'rb')
	#data = pickle.load(binin)
	#binin.close()

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Done')

	sys.exit(0)
