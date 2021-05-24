#!/usr/bin/env python3

import os
import sys
import pickle
import argparse
from argparse import HelpFormatter
from datetime import datetime

#additional modules

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import AgglomerativeClustering



class CustomFormat(HelpFormatter):

	'''
	Custom help format
	'''


	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()



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
	Get input and dump dendrogram and silhouette analysis to file
	'''


	parser = argparse.ArgumentParser(prog='TREADMILL', description='''Plot dendrogram and perform Silhouette analysis for agglomerative clustering''', epilog='''This program was developed by Davide Bolognini (https://github.com/davidebolo1993)''', formatter_class=CustomFormat) 

	required=parser.add_argument_group('Required I/O arguments')

	required.add_argument('-d', '--dendrogram', help='BIN file containing dendrogram map', metavar='BIN', required=True)
	required.add_argument('-s', '--similarity_matrix', help='BIN file containing the similarity matrix', metavar='BIN', required=True)
	required.add_argument('-o', '--output', help='output folder', metavar='DIR', required=True)

	additionals=parser.add_argument_group('Additional arguments')

	additionals.add_argument('--width_dendrogram', help='width of dendrogram plot [30.0]', metavar='', default=30.0, type=float)
	additionals.add_argument('--height_dendrogram', help='height of dendrogram plot [10.0]', metavar='', default=10.0, type=float)
	additionals.add_argument('--width_silhouette', help='width of silhouette plot [10.0]', metavar='', default=10.0, type=float)
	additionals.add_argument('--height_silhouette', help='height of silhouette plot [20.0]', metavar='', default=20.0, type=float)
	additionals.add_argument('--maxclust', help='maximum number of clusters for silhouette analysis [6]', metavar='', default=6, type=int)

	args = parser.parse_args()	
	
	resdict=dict()

	OUTDIR=os.path.abspath(args.output)

	if not os.path.exists(OUTDIR):

		try:

			os.makedirs(OUTDIR)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Error] Cannot create the output folder')
			sys.exit(1)

	if not os.access(OUTDIR,os.W_OK):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Missing write permissions on the output folder')
		sys.exit(1)


	OUTD=os.path.abspath(OUTDIR+ '/dendrogram.pdf')
	OUTS=os.path.abspath(OUTDIR + '/silhouette.pdf')
	OUTT=os.path.abspath(OUTDIR + '/silhouettescores.tsv')

	BIND=os.path.abspath(args.dendrogram)

	if not os.path.exists(BIND):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid dendrogram BIN file')
		sys.exit(1)


	BINS=os.path.abspath(args.similarity_matrix)

	if not os.path.exists(BINS):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Invalid similarity matrix BIN file')
		sys.exit(1)


	dendrogram=open(BIND,'rb')
	dendmodel = pickle.load(dendrogram)
	dendrogram.close()
	silhouette=np.load(BINS)


	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Extracting full dendrogram and plotting')

	linkage_matrix=LinkageMatrix(dendmodel)

	plt.figure(figsize=(args.width_dendrogram,args.height_dendrogram))
	plt.title('Hierarchical Clustering dendrogram')
	hierarchy.dendrogram(linkage_matrix)
	plt.ylabel('Cluster Dissimilarity Index')
	plt.savefig(OUTD)


	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + ']' + '[Message] Performing Silhouette analysis and plotting')

	range_n_clusters = range(2,args.maxclust+1)

	fig,axs = plt.subplots(len(range_n_clusters),2)
	fig.set_size_inches(args.width_silhouette,args.height_silhouette)

	for i,n_clusters in enumerate(range_n_clusters):

		axs[i,0].set_xlim([-1, 1])
		axs[i,0].set_ylim([0, len(silhouette) + (n_clusters + 1) * 10])
		clusterer = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed', linkage='average')
		cluster_labels = clusterer.fit_predict(silhouette)
		silhouette_avg = silhouette_score(silhouette, cluster_labels, metric='precomputed')

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message]' + ' For number of clusters = ' + str(n_clusters) + ', the average silhouette_score is : ' + str(silhouette_avg))
		resdict[str(n_clusters)] = float(silhouette_avg)

		sample_silhouette_values = silhouette_samples(silhouette, cluster_labels)
		y_lower = 10
		
		for l in range(n_clusters):

			ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == l]
			ith_cluster_silhouette_values.sort()
			size_cluster_i = ith_cluster_silhouette_values.shape[0]
			y_upper = y_lower + size_cluster_i
			color = cm.nipy_spectral(float(l) / n_clusters)
			axs[i,0].fill_betweenx(np.arange(y_lower, y_upper),0, ith_cluster_silhouette_values,facecolor=color, edgecolor=color, alpha=0.7)
			#axs[i,0].text(-0.05, y_lower + 0.5 * size_cluster_i, str(l))
			y_lower = y_upper + 10 

		axs[i,0].set_title('Silhouette plot with number of clusters = ' + str(n_clusters))
		axs[i,0].set_xlabel('Silhouette coefficient value')
		axs[i,0].set_ylabel('Clusters')
		axs[i,0].axvline(x=silhouette_avg, color="red", linestyle="--")
		axs[i,0].set_yticks([])
		axs[i,0].set_xticks([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])

		colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
		axs[i,1].scatter(silhouette[:, 0], silhouette[:, 1], marker='.', s=30, lw=0, alpha=0.7,c=colors, edgecolor='k')
		axs[i,1].set_title('Cluster plot with number of clusters = ' + str(n_clusters))
		axs[i,1].set_xlabel('Feature space (1st feature)')
		axs[i,1].set_ylabel('Feature space (2nd feature)')

	plt.tight_layout()
	plt.savefig(OUTS)

	reslist=sorted(resdict.items(), key=lambda item: item[1], reverse=True)
	tab=pd.DataFrame({'n_clusters': [x[0] for x in reslist], 'silhouette_score' : [x[1] for x in reslist]})
	tab.to_csv(OUTT, sep='\t', header=True, index=False)


if __name__ == '__main__':

	main()
