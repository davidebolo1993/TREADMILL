#!/usr/bin/python env

#python 3 standard library

import argparse
import sys
from argparse import HelpFormatter

def main():

	
	parser = argparse.ArgumentParser(prog='TREADMILL', description='''TREADMILL: Tandem REpeats AnD MethylatIon caLLing''', epilog='''This program was developed by Davide Bolognini (https://github.com/davidebolo1993)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='BASIC')

	## BASIC ##

	parser_basic = subparsers.add_parser('BASIC', help='BAm StatIstiCs. Calculate BAM statistics for on-target reads from a crispr-cas9 nanopore run')

	required = parser_basic.add_argument_group('Required I/O arguments')

	required.add_argument('-bam', '--bamfile', help='sorted and MD-tagged BAM from minimap2/NGMLR', metavar='BAM', required=True)
	required.add_argument('-bed', '--bedfile', help='on-target regions in BED format', metavar="BED", required=True)
	required.add_argument('-o', '--output', help='output JSON file', metavar='JSON', required=True)

	additional = parser_basic.add_argument_group('Additional parameters')

	additional.add_argument('-z', '--gzipped', help='output gzipped JSON file', action='store_true')

	parser_basic.set_defaults(func=run_subtool)

	#print help if no subcommand nor --help provided
	
	if len(sys.argv)==1:
    	
		parser.print_help(sys.stderr)
		sys.exit(1)

	#case-insensitive submodules
	
	if sys.argv[1].lower() == 'basic':

		sys.argv[1] = 'BASIC'

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
	
	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':


	main()