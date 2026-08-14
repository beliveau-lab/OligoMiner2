"""Tests for the command-line entry point."""

import pytest

from oligominer import __version__
from oligominer.cli.main import build_parser, main


class TestParser:

    def test_reports_the_package_version(self, capsys):
        with pytest.raises(SystemExit) as exit_info:
            build_parser().parse_args(['--version'])

        assert exit_info.value.code == 0
        assert __version__ in capsys.readouterr().out

    def test_the_program_name_is_oligominer(self):
        assert build_parser().prog == 'oligominer'

    def test_an_unknown_command_is_rejected(self):
        with pytest.raises(SystemExit) as exit_info:
            build_parser().parse_args(['not-a-command'])

        assert exit_info.value.code != 0


class TestMain:

    def test_no_command_prints_help_and_fails(self, capsys):
        status = main([])

        assert status == 1
        assert 'oligominer' in capsys.readouterr().out

    def test_a_registered_command_is_dispatched(self, monkeypatch):
        # a command module sets func on its parser; main dispatches on it and
        # returns whatever the command returns
        import argparse

        from oligominer.cli import main as cli

        def parser_with_demo():
            parser = argparse.ArgumentParser(prog='oligominer')
            subparsers = parser.add_subparsers(dest='command')
            subparsers.add_parser('demo').set_defaults(func=lambda args: 7)
            return parser

        monkeypatch.setattr(cli, 'build_parser', parser_with_demo)
        assert cli.main(['demo']) == 7
