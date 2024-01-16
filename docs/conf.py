# Project settings.
project = 'Goupil'
copyright = 'Université Clermont Auvergne, CNRS/IN2P3, LPC'
author = 'Valentin Niess'
release = '0.1.3'

highlight_language = 'python3'

# General settings.
extensions = [
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

rst_prolog = """
.. role:: python(code)
    :language: python
    :class: highlight

.. role:: cpp(code)
    :language: c++
    :class: highlight
"""

# Mappings for links to externals documentations.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None)
}

# Toctrees options.
toc_object_entries = True
toc_object_entries_show_parents = 'hide'

# HTML options.
html_theme_options = {
    'collapsiblesidebar': True,
    'root_include_title': False
}
html_sidebars = {
    '**': [
        'localtoc.html',
        'relations.html',
        'sourcelink.html',
    ]
}
html_theme = 'python_docs_theme'
html_short_title = 'Goupil 0.1.3 Documentation'
html_split_index = True
html_static_path = ['_static']

# Doctest options.
doctest_global_setup = "import goupil"
