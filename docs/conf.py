"""sphinx configuration."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oligominer import __version__

project = "oligominer"
author = "Conor Camplisson"                 
copyright = "2026, Conor Camplisson"
release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # reads the Google-style docstrings
    "myst_parser",  # docs pages are Markdown, like the README
]

# napoleon is configured for Google style only; numpy style is not used here.
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# RENDERING, not content. By default napoleon turns an `Attributes:` block into one
# py:attribute directive per entry, and each type becomes its own full-width "Type:" table —
# a table whose entire content is the word "str". use_ivar renders them as a compact field
# list on the class instead. use_rtype folds the return type into the Returns line rather
# than emitting a separate "Return type:" block.
napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = True

# "Returns: None" appears throughout this package's docstrings; it is noise in rendered output
# and says nothing a reader did not already assume.
napoleon_include_init_with_doc = False

html_theme = "pydata_sphinx_theme"
# The switcher entry for a development build is literally "dev", because that is the folder CI
# writes it to; a released build matches its own version. Without this the dropdown on the dev
# site would look for "0.2.0.dev1", find no such entry, and highlight nothing.
version_match = "dev" if "dev" in release else release

html_theme_options = {
    # the version switcher is populated by the docs workflow writing switcher.json
    "switcher": {
        "json_url": "https://oligominer.org/switcher.json",
        "version_match": version_match,
    },
    "navbar_end": ["version-switcher", "theme-switcher", "navbar-icon-links"],
    "github_url": "https://github.com/beliveau-lab/OligoMiner2",
    "use_edit_page_button": False,
}

html_title = f"{project} {release}"

exclude_patterns = ["_build"]
