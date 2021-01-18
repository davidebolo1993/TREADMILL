#!/usr/bin/python env

#python 3 standard library

import argparse
import sys
from argparse import HelpFormatter

def main():

	parser = argparse.ArgumentParser(prog='TREADMILL', description='''TREADMILL: Tandem REpeats AnD MethylatIon caLLing''', epilog='''This program was developed by Davide Bolognini (https://github.com/davidebolo1993)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='BASIC,RYDER,TRAP')

	## BASIC ##

	parser_basic = subparsers.add_parser('BASIC', help='BAm StatIstiCs. Calculate BAM statistics for on-target reads from a targeted nanopore experiment')

	required = parser_basic.add_argument_group('Required I/O arguments')

	required.add_argument('-bam', '--bamfile', help='sorted and MD-tagged BAM file', metavar='BAM', required=True)
	required.add_argument('-bed', '--bedfile', help='on-target regions in BED format', metavar='BED', required=True)
	required.add_argument('-o', '--output', help='output gzipped JSON file', metavar='JSON.GZ', required=True)

	parser_basic.set_defaults(func=run_subtool)

	## RYDER ##

	parser_reef = subparsers.add_parser('RYDER', help='gRoup bY DEcision tRee. Extract on-target reads from a targeted nanopore experiment and group them by similarity. This module (1) a re-maps the original reads to synthetic chromosomes harboring repeat expansions and (2) group reads by similarity using an edit distance-based decision tree')

	required = parser_reef.add_argument_group('Required I/O arguments')

	required.add_argument('-fa', '--fastafile', help='reference genome in FASTA format', metavar='FASTA', required=True)
	required.add_argument('-bam', '--bamfile', help='sorted BAM file', metavar='BAM', required=True)
	required.add_argument('-bed', '--bedfile', help='on-target regions in BED format', metavar='BED', required=True)
	required.add_argument('-o', '--output', help='output binary map. Output folder will be created if it does not exist', metavar='BIN', required=True)
	required.add_argument('--motif', help='known repeated motif (one for each region in the BED file given to RYDER)', nargs='+', action='append', required=True, metavar='MOTIF')


	cluster = parser_reef.add_argument_group('Agglomerative clustering parameters')
	cluster.add_argument('--affinity', help='sequence similarity percentage between grouped sequences [70.0]', type=float, default=70.0, metavar='')
	cluster.add_argument('--support', help='minimum group support (retain only clustered groups with enough members) [5]', default=5, type=int, metavar='')
	cluster.add_argument('--hierarchical_clustering', help = 'perform agglomerative hierarchical clustering instead of using the greedy approach. One between --threshold and --clusters must be specified', action='store_true')
	cluster.add_argument('--treshold', help = 'cut dendogram at given treshold [None]', default=None, metavar='')
	cluster.add_argument('--clusters', help = 'output specified number of clusters [None]', default=None, metavar='')
	

	additional = parser_reef.add_argument_group('Additional parameters')

	additional.add_argument('--maxsize', help='maximum number (approximate) of repeated motifs in the (synthetic) reference sequences [500]', type=int, default=500, metavar='')
	additional.add_argument('--flanking', help='number of bases flanking repeats in the (synthetic) reference sequences [1000]', type=int, metavar='', default=1000)
	additional.add_argument('--similarity', help='sequence similarity percentage between generated (synthetic) reference sequences [85.0]', type=float, metavar='', default=85.0)
	additional.add_argument('--threads', help='number of threads to use for the re-alignment step [1]', type=int, metavar='', default=1)
	additional.add_argument('--store', help='store the synthetic chromosomes used for the re-alignment step in FASTA file and the re-aligned BAM in the same folder used for the BIN file', action='store_true')
	
	parser_reef.set_defaults(func=run_subtool)

	## TRAP ##

	parser_trap = subparsers.add_parser('TRAP', help='Tandem RepeAts Profiler. Identify and genotype tandem repeats from clusters of reads created with RYDER')

	required = parser_trap.add_argument_group('Required I/O arguments')

	required.add_argument('-i', '--input', help='input binary map from RYDER', metavar='BIN', required=True)
	required.add_argument('-o', '--output', help='output directory', metavar='DIR', required=True)

	algorithm = parser_trap.add_argument_group('Consensus sequence computation')

	algorithm.add_argument('-m', '--match', help='match reward for consensus computation [5]', metavar='', default=5, type=int)
	algorithm.add_argument('-x', '--mismatch', help='mismatch penalty for consensus computation [-4]', metavar='', default=-4, type=int)
	algorithm.add_argument('-go', '--gapopen', help='gap opening penalty for consensus computation [-8]', metavar='', default=-8, type=int)
	algorithm.add_argument('-ge', '--gapextend', help='gap extending penalty for consensus computation [-6]', metavar='', default=-6, type=int)
	algorithm.add_argument('--reference_weight', help='reference weight for consensus computation. This weights the reference with respect to the number of noisy reads for which consensus is computed [0.0]', metavar='', default=0.0, type=float)

	repeat=parser_trap.add_argument_group('Repeat profiling parameters')

	repeat.add_argument('--substitution', help='substitution cost for weighted edit distance calculation [1.0]', metavar='', default=1.0, type=float)
	repeat.add_argument('--deletion', help='deletion cost for weighted edit distance calculation [1.0]', metavar='', default=1.0, type=float)
	repeat.add_argument('--insertion', help='insertion cost for weighted edit distance calculation [1.0]', metavar='', default=1.0, type=float)
	repeat.add_argument('--maxedit', help='maximum edit distance score for approximate string matching. If the edit distance between the repeated motif and a substring is greater than this value, then this is ignored. Otherwise the substring is considered either a repeat interruption or an approximate repeat [1.0]', metavar='', default=1.0, type=float)

	additional = parser_trap.add_argument_group('Additional parameters')

	additional.add_argument('--similarity', help='sequence similarity percentage (discriminate between reference and alternative alleles) [80.0]', required=False, default=80.0, type=float, metavar='')
	additional.add_argument('--subgroups', help='if multiple alleles are present, output their repeat content and frequency', action='store_true')

	parser_trap.set_defaults(func=run_subtool)

	#print help if no subcommand nor --help provided

	print(r"""

	 ______ ______  ______  ______  _____   __    __  __  __      __        
	/\__  _/\  == \/\  ___\/\  __ \/\  __-./\ "-./  \/\ \/\ \    /\ \       
	\/_/\ \\ \  __<\ \  __\\ \  __ \ \ \/\ \ \ \-./\ \ \ \ \ \___\ \ \____  
	   \ \_\\ \_\ \_\ \_____\ \_\ \_\ \____-\ \_\ \ \_\ \_\ \_____\ \_____\ 
	    \/_/ \/_/ /_/\/_____/\/_/\/_/\/____/ \/_/  \/_/\/_/\/_____/\/_____/ v1.0
	                                                                        
	""")
	
	
	if len(sys.argv)==1:
    	
		parser.print_help(sys.stderr)
		sys.exit(1)

	#case-insensitive submodules
	
	if sys.argv[1].lower() == 'basic':

		sys.argv[1] = 'BASIC'

	elif sys.argv[1].lower() == 'ryder':

		sys.argv[1] = 'RYDER'

	elif sys.argv[1].lower() == 'trap':

		sys.argv[1] = 'TRAP'

	args = parser.parse_args()
	args.func(parser, args)


class CustomFormat(HelpFormatter):


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



def run_subtool(parser, args):


	if args.command == 'BASIC': #BAm StatIstiCs

		from .BASIC import BASIC as submodule

	elif args.command == 'RYDER': #gRoup bY DEcision tRee

		from .RYDER import RYDER as submodule

	elif args.command == 'TRAP': #Tandem RepeAts Profiler 

		from .TRAP import TRAP as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':


	main()
