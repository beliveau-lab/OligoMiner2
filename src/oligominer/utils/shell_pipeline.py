"""
# Shell Pipeline

Utilities for running shell commands from Python. Use run_cmd() for single
commands, or ShellPipeline for chaining multiple commands together with
stdout-to-stdin piping between steps. Both support feeding input data from
RAM and writing final output to a file.
"""

import subprocess
import sys

from .exceptions import ExternalCommandFailed, PipelineStateError, MissingOutputFile


def run_cmd(cmd, input_data=None, output_file=None, binary=False, verbose=False):
    """
    Run a single shell command.

    Convenience wrapper around ShellPipeline for the common case of running
    one command without chaining.

    Args:
        cmd (list): command and its arguments as a list.
        input_data (str or bytes, optional): data to pipe into stdin.
        output_file (str, optional): if provided, write output to this file.
        binary (bool): if True, operate in binary mode; otherwise, text mode.
        verbose (bool): if True, print stdout and stderr to the terminal.

    Returns:
        result (str, bytes, or None): command output, or None if output_file
            is provided.

    Raises:
        ExternalCommandFailed: if the command fails.
    """
    result = ShellPipeline(binary=binary).add(cmd).run(
        input_data=input_data, output_file=output_file, verbose=verbose
    )

    # success
    return result

class ShellPipeline:
    def __init__(self, binary=False):
        """
        Initialize the pipeline.

        Args:
            binary (bool): if True, operate in binary mode; otherwise, text mode.
        """
        self.commands = []
        self.binary = binary

    def add(self, cmd):
        """
        Add a command to the pipeline.

        Args:
            cmd (list): command and its arguments as a list.

        Returns:
            self (ShellPipeline): for method chaining.
        """
        self.commands.append(cmd)
        return self

    def add_multi(self, *cmds):
        """
        Add multiple commands to the pipeline.

        Args:
            *cmds (list): commands and their arguments as lists.

        Returns:
            self (ShellPipeline): for method chaining.
        """
        for cmd in cmds:
            self.add(cmd)
        return self

    def run(self, input_data=None, output_file=None, verbose=False):
        """
        Run the entire pipeline.

        Args:
            input_data (str or bytes, optional): data to pass to the first command.
            output_file (str, optional): if provided, final output will be written to this file.
            verbose (bool): if True, print stdout and stderr of each command to the terminal.

        Returns:
            result (str, bytes, or None): the output of the final command, or None if
                output_file is provided.

        Raises:
            PipelineStateError: if no commands have been added to the pipeline.
            ExternalCommandFailed: if any command in the pipeline fails.
        """
        if not self.commands:
            raise PipelineStateError("No commands in the pipeline.")

        text_mode = not self.binary
        current_stdout = self._run_command(self.commands[0], input_data, text_mode, verbose)

        for cmd in self.commands[1:]:
            current_stdout = self._run_command(cmd, current_stdout, text_mode, verbose)

        if output_file:
            self._write_output(output_file, current_stdout)
            return None

        # success
        return current_stdout

    def _run_command(self, cmd, input_data, text_mode, verbose):
        """
        Run a single command in the pipeline.

        Args:
            cmd (list): command and its arguments as a list.
            input_data (str or bytes, optional): data to pass to the command.
            text_mode (bool): if True, operate in text mode; otherwise, binary mode.
            verbose (bool): if True, print stdout and stderr of the command to the terminal.

        Returns:
            current_stdout (str or bytes): the output of the command.

        Raises:
            ExternalCommandFailed: if the command fails to start or returns
                a non-zero exit code.
        """
        try:
            with subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE if input_data is not None else None,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=text_mode
            ) as p:
                current_stdout, err = p.communicate(input=input_data)
        except Exception as e:
            raise ExternalCommandFailed(cmd, returncode=-1, stderr=str(e))

        if verbose:
            if current_stdout:
                sys.stdout.write(current_stdout)
            if err:
                sys.stderr.write(err)
            sys.stdout.flush()
            sys.stderr.flush()

        if p.returncode != 0:
            raise ExternalCommandFailed(cmd, returncode=p.returncode, stderr=err.strip() if err else None)

        # success
        return current_stdout

    def _write_output(self, output_file, data):
        """
        Write the final output to a file.

        Args:
            output_file (str): the file to write the output to.
            data (str or bytes): the data to write.

        Raises:
            MissingOutputFile: if writing to the file fails.
        """
        mode = 'wb' if self.binary else 'w'
        try:
            with open(output_file, mode) as f:
                f.write(data)
        except Exception as e:
            raise MissingOutputFile(output_file) from e
