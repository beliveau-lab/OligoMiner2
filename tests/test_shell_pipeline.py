"""
Tests for running external commands.

Every external tool in the package goes through here, so a failure that does
not raise becomes a silently wrong result downstream. The failure paths are
tested as carefully as the success ones.
"""

import pytest

from oligominer.utils.exceptions import (
    ExternalCommandFailed, PipelineStateError,
)
from oligominer.utils.shell_pipeline import ShellPipeline, run_cmd


class TestRunCmd:

    def test_returns_the_command_output(self):
        assert run_cmd(['echo', 'hello']).strip() == 'hello'

    def test_input_is_piped_to_stdin(self):
        assert run_cmd(['cat'], input_data='piped').strip() == 'piped'

    def test_output_goes_to_a_file_when_asked(self, tmp_path):
        path = tmp_path / 'out.txt'

        assert run_cmd(['echo', 'written'], output_file=str(path)) is None
        assert path.read_text().strip() == 'written'

    def test_a_failing_command_raises(self):
        with pytest.raises(ExternalCommandFailed):
            run_cmd(['false'])

    def test_a_missing_executable_raises(self):
        with pytest.raises(ExternalCommandFailed):
            run_cmd(['this-command-does-not-exist'])

    def test_binary_mode_returns_bytes(self):
        assert isinstance(run_cmd(['echo', 'bytes'], binary=True), bytes)


class TestShellPipeline:

    def test_stages_are_chained(self):
        result = (ShellPipeline()
                  .add(['printf', 'b\na\nc\n'])
                  .add(['sort'])
                  .run())

        assert result.split() == ['a', 'b', 'c']

    def test_add_multi_chains_in_order(self):
        result = (ShellPipeline()
                  .add_multi(['printf', 'one\ntwo\n'], ['wc', '-l'])
                  .run())

        assert result.strip() == '2'

    def test_add_returns_self_for_chaining(self):
        pipeline = ShellPipeline()
        assert pipeline.add(['true']) is pipeline

    def test_an_empty_pipeline_raises(self):
        with pytest.raises(PipelineStateError):
            ShellPipeline().run()

    def test_a_failure_in_a_later_stage_raises(self):
        # a stage that fails after data has already flowed must not be
        # mistaken for an empty result
        with pytest.raises(ExternalCommandFailed):
            ShellPipeline().add(['echo', 'data']).add(['false']).run()

    def test_a_failure_in_the_first_stage_raises(self):
        with pytest.raises(ExternalCommandFailed):
            ShellPipeline().add(['false']).add(['cat']).run()

    def test_input_reaches_the_first_stage(self):
        result = (ShellPipeline()
                  .add(['cat'])
                  .add(['tr', 'a-z', 'A-Z'])
                  .run(input_data='shout'))

        assert result.strip() == 'SHOUT'

    def test_output_file_returns_none_and_writes(self, tmp_path):
        path = tmp_path / 'piped.txt'
        result = (ShellPipeline()
                  .add(['printf', 'x\ny\n'])
                  .add(['wc', '-l'])
                  .run(output_file=str(path)))

        assert result is None
        assert path.read_text().strip() == '2'
