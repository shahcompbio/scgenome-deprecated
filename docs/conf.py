# Configuration file for the Sphinx documentation builder.

import os
import sys
from pathlib import Path

from packaging.version import parse as parse_version

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent)]

import scgenome  # noqa

on_rtd = os.environ.get('READTHEDOCS') == 'True'

# -- General configuration ------------------------------------------------


nitpicky = True  # Warn about broken links. This is here for a reason: Do not change.
needs_sphinx = '2.0'  # Nicer param docs
suppress_warnings = [
    'ref.citation',
    'myst.header',  # https://github.com/executablebooks/MyST-Parser/issues/262
]

project = 'scgenome'
copyright = '2018, McPherson'
author = 'Andrew McPherson'

release = '0.0'
version = '0.0.7'



# Bumping the version updates all docs, so don't do that
if parse_version(version).is_devrelease:
    parsed = parse_version(version)
    version = f"{parsed.major}.{parsed.minor}.{parsed.micro}.dev"

# default settings
templates_path = ['_templates']
master_doc = 'index'
default_role = 'literal'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
source_suffix = [".rst", ".md"]

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'matplotlib.sphinxext.plot_directive',
    'nbsphinx',
    'sphinx_gallery.load_style',
    'sphinx_autodoc_typehints',  # needs to be after napoleon
    'scanpydoc.rtd_github_links',
]


# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False

typehints_defaults = 'braces'

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}

# -- Options for HTML output ----------------------------------------------

html_theme = 'scanpydoc'
html_context = {
    'display_github': True,
    "github_user": "shahcompbio",
    "github_repo": "scgenome",
    "github_version": "master",
    "conf_py_path": "/docs/",
}


def setup(app):
    app.warningiserror = on_rtd


# -- Options for other output formats ------------------------------------------

htmlhelp_basename = f'{project}doc'
doc_title = f'{project} Documentation'
latex_documents = [(master_doc, f'{project}.tex', doc_title, author, 'manual')]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (
        master_doc,
        project,
        doc_title,
        author,
        project,
        'One line description of project.',
        'Miscellaneous',
    )
]


# Options for plot examples
nitpick_ignore = [
    ('py:class', 'type'),
    ('py:class', 'AnnData'),
    ('py:class', 'Axes'),
    ('py:class', 'Figure'),
    ('py:class', 'optional'),
    ('py:class', 'Bio.Phylo.BaseTree.Tree'),
    ('py:class', 'matplotlib.colors.ListedColormap'),
    ('py:class', 'anndata.AnnData'),
    ('py:class', 'ad.AnnData'),
    ('py:class', 'anndata._core.anndata.AnnData'),
    ('py:class', 'DataFrame'),
    ('py:class', 'Sequence'),
    ('py:class', 'PyRanges'),
    ('py:class', 'pyranges.PyRanges'),
    ('py:class', 'pyrange.PyRanges'),
    ('py:class', 'pyranges.pyranges.PyRanges'),
    ('py:class', 'ndarray'),
    ('py:class', 'matplotlib.axes.Axes'),
    ('py:class', 'matplotlib.figure.Figure'),
    # Currently undocumented: https://github.com/mwaskom/seaborn/issues/1810
    ('py:class', 'seaborn.ClusterGrid'),
    ('py:class', 'numpy.random.mtrand.RandomState'),
    # Will work once scipy 1.8 is released
    ('py:class', 'scipy.sparse.base.spmatrix'),
    ('py:class', 'scipy.sparse.csr.csr_matrix'),
]

plot_include_source = True
plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_working_directory = HERE.parent  # Project root
