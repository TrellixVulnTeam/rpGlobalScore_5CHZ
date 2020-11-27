.. rpGlobalScore documentation master file, created by
   sphinx-quickstart on Fri Nov 27 11:03:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to rpGlobalScore's documentation!
=========================================

Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
############

.. _rpBase: https://github.com/Galaxy-SynBioCAD/rpBase

Welcome to the documentation of the rpGlobalScore project. This project combines the different criteria calculated by different tools of the same collection of tools, including: Thermodynamics, FBA, Reaction Rule score, and length of the heterologous pathway. These criterias are combined using a weighted sum, where the weights have been optimised to 81 experimentally implemented heterologous pathways compared to the results of the pipeline. These weights can be changed by the user, but must be between 0.0 and 1.0.

Usage
#####

First build the rpBase_ docker before building the local docker using the following command:

.. code-block:: bash

   docker build -t brsynth/rpglobalscore-standalone -f Dockerfile .

To run the docker in the command line, please use the run.py script. An example input could be:

.. code-block:: bash

   python run.py -input /path/to/file_rpsbml.tar.gz -input_format tar -output /path/to/outfile.tar.gz

API
###

.. toctree::
   :maxdepth: 2
   :caption: Contents:



