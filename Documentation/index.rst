.. title:: Home

.. meta::
   :description: Overview over the Medical Image Registration ToolKit (MIRTK)
   :keywords:    image processing, image registration, IRTK, MIRTK

==========================================
Medical Image Registration ToolKit (MIRTK)
==========================================

The MIRTK is a research-focused image processing toolkit, developed at the
BioMedIA_ research group. It provides a collection of libraries and command-line
tools to assist in processing and analyzing imaging data.
The main application of the MIRTK is in adult and neonatal brain MR image registration
as well as the reconstruction of cortical surface meshes. The modular project
organization of the MIRTK enables the installation of selected modules.

In the event you found the MIRTK useful, please consider giving appropriate credit
to the software with a citation of the research article(s) describing the implemented
algorithm(s). See the :doc:`list of publications <publications>` for suitable references.


Modules
=======

.. include:: modules/_overview.rst


Commands
========

.. include:: commands/_overview.rst


Background
==========

Parts of the Common, Numerics, Image, Transformation, and Registration modules and
command-line tools of the MIRTK originated from the `IRTK`_ written by `Daniel Rueckert`_
and `Julia Schnabel`_. All of the transformation and registration code of the IRTK was
rewritten from scratch by `Andreas Schuh`_ during his PhD studies, with a new modular and
extended registration framework. Additional modules and commands for the reconstruction
and inflation of cortical surface meshes, the registration of surface meshes, and the
harmonic mapping of brain volumes were subsequently added to the MIRTK.

.. _BioMedIA: https://biomedia.doc.ic.ac.uk/
.. _Daniel Rueckert: http://www.imperial.ac.uk/people/d.rueckert
.. _Julia Schnabel: http://www.imagingcdt.com/about/cdt-management/professor-julia-schnabel
.. _Andreas Schuh: http://www.andreasschuh.com/
.. _IRTK: https://github.com/BioMedIA/IRTK
