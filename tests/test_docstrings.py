"""Tests that the public surface is documented.

The package documentation is generated from docstrings, so a public function
without one, or without its Args and Returns, is a hole in the docs rather than a
style lapse.

The vendored numerical kernels under l4t/ and bilstm_arch.py are excluded: they
are copied verbatim from the studies that fitted the shipping models, and
reformatting them would risk changing the model's predictions for no
documentation gain.
"""

import ast
import pathlib

import pytest

import oligominer

PACKAGE_ROOT = pathlib.Path(oligominer.__file__).parent

# copied verbatim so the shipped models' predictions do not move
VENDORED = ('l4t', 'bilstm_arch.py')


def public_definitions():
    """
    Yield every public function and class defined in the package.

    Returns:
        definitions (list): (path, node) pairs for each public definition.
    """
    definitions = []
    for path in sorted(PACKAGE_ROOT.rglob('*.py')):
        if any(part in VENDORED for part in path.parts) or path.name in VENDORED:
            continue
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                if not node.name.startswith('_'):
                    definitions.append((path, node))

    # success
    return definitions


DEFINITIONS = public_definitions()
IDS = [f'{p.relative_to(PACKAGE_ROOT)}::{n.name}' for p, n in DEFINITIONS]


def test_the_package_has_public_definitions_to_check():
    assert len(DEFINITIONS) > 100


@pytest.mark.parametrize('path,node', DEFINITIONS, ids=IDS)
def test_has_a_docstring(path, node):
    assert (ast.get_docstring(node) or '').strip(), (
        f'{node.name} in {path.name} has no docstring, so it will be undocumented')


@pytest.mark.parametrize('path,node', DEFINITIONS, ids=IDS)
def test_documents_its_arguments(path, node):
    if not isinstance(node, ast.FunctionDef):
        return
    doc = ast.get_docstring(node) or ''
    args = [a.arg for a in node.args.args if a.arg not in ('self', 'cls')]
    if args and doc.strip():
        assert 'Args:' in doc, (
            f'{node.name} in {path.name} takes {args} but documents no Args section')


@pytest.mark.parametrize('path,node', DEFINITIONS, ids=IDS)
def test_documents_what_it_returns(path, node):
    if not isinstance(node, ast.FunctionDef):
        return
    doc = ast.get_docstring(node) or ''
    returns_a_value = any(isinstance(n, ast.Return) and n.value is not None
                          for n in ast.walk(node))
    if returns_a_value and doc.strip():
        assert 'Returns:' in doc, (
            f'{node.name} in {path.name} returns a value but documents no '
            f'Returns section')
