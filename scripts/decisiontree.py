#!/usr/bin/env python3

#python 3 standard library

from itertools import combinations

#additional modules

import editdistance


def similarity(worda,wordb):

	'''
	Return the edit-distance based similarity score between 2 sequences
	'''

	return 100-100*editdistance.eval(worda,wordb)/(len(worda)+len(wordb))


def decisiontree(wordsdict, nameslist,mingroupsize=1, treshold=85.0):

	'''
	Group strings in list by similarity (edit distance score)
	'''
	
	#wordsdict=dict(zip(rnames, sequences))

	paired ={c:{c} for c in wordsdict.values()}

	for worda,wordb in combinations(wordsdict.values(),2):

		if similarity(worda,wordb) < treshold: 

				continue

		else:

			paired[worda].add(wordb)
			paired[wordb].add(worda)


	decision = list()
	ungrouped = set(wordsdict.values())
	
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


	#retrive read names

	winners = list()

	for groups in decision:

		groupkeys = set()

		for elements in groups: #retrieve all the original keys

			keys=set([k for k,v in wordsdict.items() if v == elements])

			groupkeys.update(keys)

		winners.append(groupkeys)


	return decision,winners
