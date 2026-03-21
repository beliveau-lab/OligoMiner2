"""
# OligoMiner Command-Line Interface (CLI)

This module sets up the command-line interface for OligoMiner, allowing 
users to run various commands related to oligonucleotide probe design 
directly from the terminal.

Commands are organized into subcommands, each handled by its own module. 
This modular design makes it easy to add new functionality in the future.
"""


import argparse

from .. import __version__
from . import test_cli # import each command module

# configure main CLI help text
HELP_TEXT = f'''
  ____  _ _             __  __ _               _____ _____ 
 / __ \| (_)           |  \/  (_)             |_   _|_   _|
| |  | | |_  __ _  ___ | \  / |_ _ __   ___ _ __| |   | |  
| |  | | | |/ _` |/ _ \| |\/| | | '_ \ / _ \ '__| |   | |  
| |__| | | | (_| | (_) | |  | | | | | |  __/ | _| |_ _| |_ 
 \____/|_|_|\__, |\___/|_|  |_|_|_| |_|\___|_||_____|_____|
             __/ |                                         
            |___/                                          

Version:   {__version__}
Docs:      https://oligominer.org/docs/{__version__}/
Code:      https://github.com/beliveau-lab/OligoMiner2
'''


# HELP_TEXT = f'''
# \b
#   ____  _ _             __  __ _               _____ _____ 
#  / __ \| (_)           |  \/  (_)             |_   _|_   _|
# | |  | | |_  __ _  ___ | \  / |_ _ __   ___ _ __| |   | |  
# | |  | | | |/ _` |/ _ \| |\/| | | '_ \ / _ \ '__| |   | |  
# | |__| | | | (_| | (_) | |  | | | | | |  __/ | _| |_ _| |_ 
#  \____/|_|_|\__, |\___/|_|  |_|_|_| |_|\___|_||_____|_____|
#              __/ |                                         
#             |___/                                          

# \b
# Version:   {__version__}
# Docs:      https://oligominer.org/docs/{__version__}/
# Code:      https://github.com/beliveau-lab/OligoMiner2
# '''

def build_parser():
    parser = argparse.ArgumentParser(prog="oligominer", description=HELP_TEXT, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="show version and exit",
    )

    subparsers = parser.add_subparsers(dest="command")

    # Each command module registers itself with the subparsers
    # test_cli.register(subparsers)
    # align.register(subparsers) # TODO

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if hasattr(args, "func"):
        return args.func(args)

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
