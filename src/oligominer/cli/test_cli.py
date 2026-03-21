"""
# OligoMiner CLI: `test_cli`

This module defines the `test_cli` command for the OligoMiner command-line interface (CLI).
"""


def register(subparsers):
    """
    Register the `test_cli` command with the given subparsers.

    Args:
        subparsers (argparse._SubParsersAction): The subparsers object to register the command with.

    Returns:
        None
    """
    # define parser for the CLI command
    parser = subparsers.add_parser(
        "test_cli",
        help="run a simple test command",
    )

    # set the function to be called when this command is invoked
    parser.set_defaults(func=run)


def run(args):
    """
    Run the `test_cli` command.

    Args:
        args (argparse.Namespace): The command-line arguments.

    Returns:
        exit_code (int): The exit code.
    """
    # simple test action
    print("Hello from oligominer test_cli!")

    # success
    exit_code = 0
    return exit_code
