#!/usr/bin/env python3

import os
import sys
import pickle
from datetime import datetime

#additional modules

import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy


def LinkageMatrix(model):

	'''
	Create scipy linkage matrix from sklearn model
	'''

	counts = np.zeros(model.children_.shape[0])
	n_samples = len(model.labels_)

	for i, merge in enumerate(model.children_):

		current_count = 0

		for child_idx in merge:

			if child_idx < n_samples:

				current_count += 1  # leaf node

			else:

				current_count += counts[child_idx - n_samples]

		counts[i] = current_count

	linkage_matrix = np.column_stack([model.children_, model.distances_,counts]).astype(float)

	return linkage_matrix



def main():


	'''
	Get input and dump dendogram to file
	'''
	
	BIN=os.path.abspath(sys.argv[1])
	OUT=os.path.abspath(sys.argv[2])
	OUTDIR=os.path.dirname(os.path.abspath(OUT))

	try:

		width=int(sys.argv[3])

	except:

		width=30


	try:

		height=int(sys.argv[4])

	except:

		height=10

	if not os.path.exists(BIN):

		now=datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid BIN file')
		sys.exit(1)

	if not os.path.exists(OUTDIR):

		try:

			os.makedirs(OUTDIR)

		except:

			now=datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Error] Cannot create the output folder')
			sys.exit(1)

	if not os.access(os.path.dirname(OUT),os.W_OK):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Missing write permissions on the output folder')
		sys.exit(1)


	binin=open(BIN,'rb')
	dendmodel = pickle.load(binin)
	binin.close()

	linkage_matrix=LinkageMatrix(dendmodel)

	plt.figure(figsize=(width,height))
	plt.title('Hierarchical Clustering Dendrogram')
	hierarchy.dendrogram(linkage_matrix)
	plt.savefig(OUT)


if __name__ == '__main__':

	main()
