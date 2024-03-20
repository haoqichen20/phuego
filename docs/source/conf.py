# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'phuego'
copyright = '2024, Girolamo Giudice, Haoqi Chen, Thodoris Koutsandreas, Evangelia Petsalaki'
author = 'Girolamo Giudice, Haoqi Chen, Thodoris Koutsandreas, Evangelia Petsalaki'
release = '1.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# The name of an image file (relative to this directory) to place at the top 
# of the sidebar.
html_logo = "_static/logo.png"
