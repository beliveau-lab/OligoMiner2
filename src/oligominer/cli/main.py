"""
# OligoMiner Command-Line Interface (CLI)

This module sets up the command-line interface for OligoMiner, allowing 
users to run various commands related to oligonucleotide probe design 
directly from the terminal.

Commands are organized into subcommands, each handled by its own module. A
command module exposes `register(subparsers)`, which adds its parser and sets
`func` to the callable that runs it; `main` dispatches on that attribute.
"""


import argparse

from .. import __version__

# configure main CLI help text
HELP_TEXT = rf'''
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


def build_parser():
    """
    Build the command-line argument parser.

    Each command module registers its own subparser, so adding a command does not
    require editing this function.

    Returns:
        parser (argparse.ArgumentParser): the configured parser.
    """
    parser = argparse.ArgumentParser(prog="oligominer", description=HELP_TEXT, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="show version and exit",
    )

    parser.add_subparsers(dest="command")

    return parser


def main(argv=None):
    """
    Run the command line interface.

    Args:
        argv (list, optional): arguments to parse. Defaults to sys.argv.

    Returns:
        status (int): the process exit status. 1 when no command was given, in
            which case the help text is printed.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if hasattr(args, "func"):
        return args.func(args)

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
