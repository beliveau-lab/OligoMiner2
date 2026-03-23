
from pathlib import Path
import importlib
import pkgutil

import mkdocs_gen_files

# from oligominer import alignment, config, core, duplex_stability, kmer_counting, seq_io, thermodynamics, utils

from oligominer import bioinformatics, probe_design, specificity, thermodynamics, utils

# path to this script
SCRIPT_PATH = 'docs/website/_scripts/gen_api_pages.py'

# specify which packages to document
PACKAGES = [
    bioinformatics,
    probe_design,
    specificity,
    thermodynamics,
    utils,
]

def main():
    """Main function to generate API documentation pages."""

    # generate API documentation pages for each package
    nav = mkdocs_gen_files.Nav()
    for package in PACKAGES:
        nav = build_package_docs(package=package, nav=nav)

    # write the nav file that literate-nav will read
    with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:
        # tell Material (mkdocs plugin): don't index this page
        print('---\nsearch:\n  exclude: true\n---\n', file=nav_file)
        nav_file.writelines(nav.build_literate_nav())






def build_package_docs(package, nav):
    """Generate API documentation pages for a given package.
    
    Args:
        package (module): The package to document.
        nav (mkdocs_gen_files.Nav): The navigation object to update.
    
    Returns:
        mkdocs_gen_files.Nav: The updated navigation object.
    """
    # dynamically discover modules in the package
    modules = discover_modules(package)

    # e.g. oligominer.utils -> utils
    rel_pkg = package.__name__.replace("oligominer.", "")

    # # 1) write package index.md that points to __init__.py docstring
    # pkg_index_path = Path("reference", rel_pkg, "index.md")
    # with mkdocs_gen_files.open(pkg_index_path, "w") as fd:
    #     print(f"# {package.__name__}\n", file=fd)
    #     print(f"::: {package.__name__}", file=fd)

    # # 2) add the package itself to the nav (so folder title is clickable)
    # nav[rel_pkg] = pkg_index_path.relative_to("reference").as_posix()

    # 1) write package index.md that points to __init__.py docstring
    pkg_index_path = Path("reference", rel_pkg, "index.md")
    with mkdocs_gen_files.open(pkg_index_path, "w") as fd:
        print(f"# {package.__name__}\n", file=fd)
        print(f"::: {package.__name__}", file=fd)

    # 2) register the index as the FIRST CHILD of the section
    #    (Material promotes this to a clickable section header)
    nav[rel_pkg, "index"] = pkg_index_path.relative_to("reference").as_posix()

    # 3) now write each module under that package
    for subpackage, module_names in modules.items():
        # "oligominer.utils" -> "utils", "oligominer.specificity.alignment" -> "specificity.alignment"
        rel_subpackage = subpackage.replace("oligominer.", "")

        # split dotted name into path/nav components
        # e.g. "specificity.alignment" -> ("specificity", "alignment")
        parts = tuple(rel_subpackage.split("."))

        for module in module_names:
            full_module_name = f"{subpackage}.{module}"

            # docs/reference/specificity/alignment/bowtie_align.md
            doc_path = Path("reference", *parts, f"{module}.md")

            with mkdocs_gen_files.open(doc_path, "w") as fd:
                print(f"# {full_module_name}\n", file=fd)
                print(f"::: {full_module_name}", file=fd)

            # build nav like: API Reference > specificity > alignment > bowtie_align
            nav[(*parts, module)] = doc_path.relative_to("reference").as_posix()

            # optional: point "edit this page" to the script
            mkdocs_gen_files.set_edit_path(doc_path, SCRIPT_PATH)

    # success
    return nav


def discover_modules(package):
    """Recursively discover all modules in a package and its subpackages.

    Args:
        package (module): the top-level package to scan.

    Returns:
        discovered (dict): mapping of dotted package name to list of module
            names, e.g. {"oligominer.utils": ["file_paths", "seq_utils"]}.
    """
    discovered = {}
    prefix = package.__name__ + "."

    for _, modname, ispkg in pkgutil.iter_modules(package.__path__, prefix):
        if ispkg:
            # recurse into subpackage
            subpkg = importlib.import_module(modname)
            sub_discovered = discover_modules(subpkg)
            discovered.update(sub_discovered)
        else:
            # leaf module — group under its parent package
            parent = prefix[:-1]
            if parent not in discovered:
                discovered[parent] = []
            discovered[parent].append(modname.split(".")[-1])

    return discovered

# run the main function
main()
