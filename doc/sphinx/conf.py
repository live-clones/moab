# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import subprocess, os, sys

def configureDoxyfile(input_dir, output_dir):
    with open('Doxyfile.in', 'r') as file :
        filedata = file.read()

    filedata = filedata.replace('@DOXYGEN_INPUT_DIR@', input_dir)
    filedata = filedata.replace('@DOXYGEN_OUTPUT_DIR@', output_dir)

    with open('Doxyfile', 'w') as file:
        file.write(filedata)

# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

breathe_projects = {}

if read_the_docs_build:
    input_dir = 'doc'
    output_dir = 'doc/sphinx'
    #configureDoxyfile(input_dir, output_dir)
    #subprocess.call('doxygen', shell=True)
    breathe_projects['MOAB'] = output_dir + '/xml'

# -- Project information -----------------------------------------------------

project = 'MOAB'
copyright = '2021, Vijay Mahadevan, Iulian Grindeanu, Rajeev Jain, Patrick Shriwise, Paul Wilson'
author = 'Vijay Mahadevan, Iulian Grindeanu, Rajeev Jain, Patrick Shriwise, Paul Wilson'

# The full version, including alpha/beta/rc tags
release = '5.2.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx_sitemap',
    'sphinx.ext.inheritance_diagram',
    'breathe',
    'exhale'
    #'myst-parser'
]

# Breathe Configuration
breathe_default_project = "MOAB"

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./src",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "MOAB Library API",
    "doxygenStripFromPath":  "..",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": False,
    #"exhaleDoxygenStdin":    "INPUT = ../src ../itaps ../tools ../pymoab ../examples"
    "exhaleDoxygenStdin":    "INPUT = ../../src"
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["*mesquite*"]

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    
    'logo_only': False,

    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}
# html_logo = ''
# github_url = ''
# html_baseurl = ''

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Breathe configuration -------------------------------------------------

#breathe_projects = {
#	"MOAB": "xml/"
#}
#breathe_default_project = "C++ Sphinx Doxygen Breathe"
#breathe_default_members = ('members', 'undoc-members')

