"""
# Custom OligoMiner2 exceptions

OligoMiner2 is designed for use in pipelines and larger workflows.
To facilitate error handling, we define a set of custom exceptions
that can be raised and caught as needed.
"""

DEFAULT_HELP_URL = 'https://oligominer.org'

class OligominerError(Exception):
    """Base error for OligoMiner2."""

    def __init__(self, message):
        super().__init__(message)
        self.message = message

    def __repr__(self):
        return f"{self.__class__.__name__}({self.message!r})"


class MissingDependency(OligominerError):
    """Something we need is not installed / importable.

    Attributes:
        dependency_name (str): Name of the missing dependency.
    """

    # don't print module full path in traceback
    __module__ = Exception.__module__

    def __init__(self, dependency_name, help_url=DEFAULT_HELP_URL):
        self.dependency_name = dependency_name
        self.help_url = help_url
        super().__init__(f"Missing dependency: {dependency_name}")

    def __str__(self):
        return (
            f"{self.dependency_name}\n\nA valid installation for {self.dependency_name} "
            f"could not be found.\n\nFor installation instructions, please visit:\n\n  {self.help_url}\n\nExiting..."
        )


class ExternalCommandFailed(OligominerError):
    """An external command failed to run successfully.

    Attributes:
        cmd (str): The command that was run.
        returncode (int): The return code from the command.
        stdout (str): Standard output from the command, if captured.
        stderr (str): Standard error from the command, if captured.
    """

    def __init__(self, cmd, returncode, stdout=None, stderr=None):
        self.cmd = cmd
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        super().__init__(f"{cmd!r} exited with {returncode}")

    def __str__(self):
        error_text = f"Command {self.cmd!r} failed with return code {self.returncode}."
        if self.stdout:
            error_text += f"\n\nStandard Output:\n{self.stdout}"
        if self.stderr:
            error_text += f"\n\nStandard Error:\n{self.stderr}"
        return error_text

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(cmd={self.cmd!r}, returncode={self.returncode}, "
            f"stdout={self.stdout!r}, stderr={self.stderr!r})"
        )


class MissingInputFile(OligominerError):
    """A required input file is missing.

    Attributes:
        filepath (str): Path to the missing file.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        super().__init__(f"Required input file is missing: {filepath}")

    def __repr__(self):
        return f"{self.__class__.__name__}(filepath={self.filepath})"


class MissingOutputFile(OligominerError):
    """An expected output file is missing.

    Attributes:
        filepath (str): Path to the missing file.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        super().__init__(f"Expected output file is missing: {filepath}")

    def __repr__(self):
        return f"{self.__class__.__name__}(filepath={self.filepath})"


class MissingDirectory(OligominerError):
    """A required directory is missing.

    Attributes:
        dirpath (str): Path to the missing directory.
    """

    def __init__(self, dirpath):
        self.dirpath = dirpath
        super().__init__(f"Required directory is missing: {dirpath}")

    def __repr__(self):
        return f"{self.__class__.__name__}(dirpath={self.dirpath})"


class InvalidInputError(OligominerError):
    """Invalid or conflicting user-supplied input parameters.

    Covers cases such as mutually exclusive arguments, out-of-range values,
    or missing required parameters.
    """

    def __init__(self, message):
        super().__init__(message)


class ConfigurationError(OligominerError):
    """Invalid configuration or incompatible parameter combinations.

    Raised when mining, model, or pipeline configuration values are
    internally inconsistent (e.g. exhaustive mode combined with spacing).
    """

    def __init__(self, message):
        super().__init__(message)


class PipelineStateError(OligominerError):
    """A pipeline step was called before its prerequisites were met.

    Raised when a method requires data from a prior step that has not
    been run yet (e.g. merge() before align()).
    """

    def __init__(self, message):
        super().__init__(message)


class DirectoryPathError(OligominerError):
    """A path that should be a directory is not a directory.

    Attributes:
        path (str): The path that is not a directory.
    """

    def __init__(self, path):
        self.path = path
        super().__init__(f"Expected a directory but no directory at path: {path}")

    def __repr__(self):
        return f"{self.__class__.__name__}(path={self.path})"
