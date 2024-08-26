# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'PICASO: Profiling Integrative Communities of Aggregated Single-cell Omics data'
copyright = '2024, Joppich, Hayat'
author = 'Joppich'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    "nbsphinx",
    "nbsphinx_link",
    'autoapi.extension'
]

exclude_patterns = ['_build', '**.ipynb_checkpoints']

nbsphinx_execute = 'never'
nbsphinx_prolog = """

{{env.docname}}
===============

"""

autoapi_add_toctree_entry = False
autoapi_dirs = ['./../../PICASO/']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'