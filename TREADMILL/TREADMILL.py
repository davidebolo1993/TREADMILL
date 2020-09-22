#!/usr/bin/python env

#python 3 standard library

import argparse
import sys
from argparse import HelpFormatter

def main():

	
	parser = argparse.ArgumentParser(prog='TREADMILL', description='''TREADMILL: Tandem REpeats AnD MethylatIon caLLing''', epilog='''This program was developed by Davide Bolognini (https://github.com/davidebolo1993)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='BASIC,READER,TRAP')

	## REEF ##

	#parser_reef = subparsers.add_parser('REEF', help='ReferEncE modiFier. Add synthetic chromosomes harboring repeat expansions to a given reference (prior to alignment)')

	#required = parser_reef.add_argument_group('Required I/O arguments')

	#required.add_argument('-fa', '--fastafile', help='reference genome in FASTA format', metavar='FASTA', required=True)
	#required.add_argument('-o', '--output', help='modified reference gnome in FASTA format', metavar='FASTA', required=True)
	#required.add_argument('region', help='repeat coordinates in RNAME[:STARTPOS[-ENDPOS]] format (samtools standard)', metavar='REGION', nargs=1)

	#additional = parser_reef.add_argument_group('Additional parameters')

	#additional.add_argument('--repeat', help='repeated motif in region [CGG]', type=str, default="CGG", metavar='')
	#additional.add_argument('--maxsize', help='maximum number of repeated motifs [500]', type=int, default=500, metavar='')
	#additional.add_argument('--contig', help='number of chromosomes with synthetic expansions (progressively longer) to generate [10]', type=int, default=10, metavar='')

	#parser_reef.set_defaults(func=run_subtool)

	## BASIC ##

	parser_basic = subparsers.add_parser('BASIC', help='BAm StatIstiCs. Calculate BAM statistics for on-target reads from a targeted nanopore experiment')

	required = parser_basic.add_argument_group('Required I/O arguments')

	required.add_argument('-bam', '--bamfile', help='sorted and MD-tagged BAM file from minimap2/NGMLR', metavar='BAM', required=True)
	required.add_argument('-bed', '--bedfile', help='on-target regions in BED format', metavar='BED', required=True)
	required.add_argument('-o', '--output', help='output JSON file', metavar='JSON', required=True)

	additional = parser_basic.add_argument_group('Additional parameters')

	additional.add_argument('-z', '--gzipped', help='output gzipped JSON file', action='store_true')

	parser_basic.set_defaults(func=run_subtool)

	## READER ##

	parser_reader = subparsers.add_parser('READER', help='READ ExtRactor. Extract on-target reads from a targeted nanopore experiment and group them by similarity')

	required = parser_reader.add_argument_group('Required I/O arguments')

	required.add_argument('-bam', '--bamfile', help='sorted BAM file from minimap2/NGMLR', metavar='BAM', required=True)
	required.add_argument('-bed', '--bedfile', help='repeated regions in BED format', metavar='BED', required=True)
	required.add_argument('-fa', '--fastafile', help='reference genome in FASTA format', metavar='FASTA', required=True)
	required.add_argument('-o', '--output', help='output binary map', metavar='BIN', required=True)

	additional = parser_reader.add_argument_group('Additional parameters')

	additional.add_argument('--similarity', help='sequence similarity percentage (discriminate group of reads with different repeat content) [85.0]', required=False, default=85.0, type=float, metavar='')
	additional.add_argument('--support', help='minimum group support (retain only groups with enough reads)[5]', required=False, default=5, type=int, metavar='')

	parser_reader.set_defaults(func=run_subtool)

	## TRAP ##

	parser_trap = subparsers.add_parser('TRAP', help='Tandem RepeAts Profiler. Identify and genotype tandem repeats from clusters of reads created with READER')

	required = parser_trap.add_argument_group('Required I/O arguments')

	required.add_argument('-i', '--input', help='input binary map from READER', metavar='BIN', required=True)
	required.add_argument('-o', '--output', help='output directory', metavar='DIR', required=True)


	algorithm = parser_trap.add_argument_group('Repeat profiling parameters')

	algorithm.add_argument('-m', '--match', help='match reward for consensus computation [5]', metavar='', default=5, type=int)
	algorithm.add_argument('-x', '--mismatch', help='mismatch penalty for consensus computation [-4]', metavar='', default=-4, type=int)
	algorithm.add_argument('-o', '--gapopen', help='gap opening penalty for consensus computation [-8]', metavar='', default=-8, type=int)
	algorithm.add_argument('-e', '--gapextend', help='gap extending penalty for consensus computation [-6]', metavar='', default=-6, type=int)
	algorithm.add_argument('--motif', help='known repeated motif (one for each region in the BED file given to READER)', nargs='+', action='append', required=True)
	
	additional = parser_trap.add_argument_group('Additional parameters')

	additional.add_argument('--similarity', help='sequence similarity percentage (discriminate between reference and alternative alleles)', required=False, default=95.0, type=float, metavar='')

	parser_trap.set_defaults(func=run_subtool)

	#print help if no subcommand nor --help provided
	
	if len(sys.argv)==1:
    	
		parser.print_help(sys.stderr)
		sys.exit(1)

	#case-insensitive submodules
	
	if sys.argv[1].lower() == 'basic':

		sys.argv[1] = 'BASIC'

	elif sys.argv[1].lower() == 'reader':

		sys.argv[1] = 'READER'

	elif sys.argv[1].lower() == 'trap':

		sys.argv[1] = 'TRAP'

	#elif sys.argv[1].lower() == 'reef':

		#sys.argv[1] = 'REEF'

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

	elif args.command == 'READER': #READ ExtRactor 

		from .READER import READER as submodule

	#elif args.command == 'REEF': #ReferEncE modiFier

		#from .REEF import REEF as submodule

	elif args.command == 'TRAP':

		from .TRAP import TRAP as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':


	main()
