.. meta::
    :description: Download the MIRTK software.
    :keywords:    MIRTK download, open source MIRTK, MIRTK license, MIRTK copyright


========
Download
========


.. _DownloadSources:

Source Code
===========

The source code of the Medical Image Registration ToolKit (MIRTK) is hosted on GitHub
by the `BioMedIA <https://github.com/BioMedIA/>`__ group.
The `BioMedIA/MIRTK <https://github.com/BioMedIA/MIRTK>`__ repository contains the
build configuration and source code files of the core libraries and basic MIRTK commands.
Additional packages can be downloaded from the `MIRTK <https://github.com/MIRTK>`__ group.

All releases and latest development versions can be downloaded by executing the
commands outlined below in a Terminal window. See the :doc:`changelog` for a summary
of changes in each release. When you download the source release packages (e.g.,
``.tar.gz`` or ``.zip`` archives) from the respective GitHub project pages instead,
you have to download and extract each individual MIRTK package separately and extract
the files into the respective subdirectory of the top-level MIRTK source tree.
We therefore recommend to use Git_ as follows to download the MIRTK source code
distribution.

To download the latest MIRTK source files using Git_, use the following
`git clone <https://git-scm.com/docs/git-clone>`__ command:

.. code-block:: bash

    git clone --branch master --depth 1 -- https://github.com/BioMedIA/MIRTK.git
    MIRTK_SOURCE_DIR="$PWD/MIRTK"

where the ``--branch master --depth 1`` options request a shallow clone of the current
stable version of the source files excluding previous revisions in order to save bandwidth,
download time, and disk space. The ``--depth`` option may be omitted if you want to
see the development history and contribute to the MIRTK. To download the latest
development version, use the ``--branch develop`` option. Note that branches can be
switched also after the initial clone command using the
`git checkout <https://git-scm.com/docs/git-checkout>`__ command. To update an existing
copy of the BioMedIA/MIRTK repository, use the `git pull <https://git-scm.com/docs/git-pull>`__
command. Note that submodules also have to be updated in this case using the commands below.

Before running any of the following Git commands, you have to change to the root directory
of the downloaded MIRTK repository:

.. code-block:: bash

    cd $MIRTK_SOURCE_DIR

.. note::

   To download the entire MIRTK distribution including
   third-party library modules, and external packages, run the command::

       git submodule update --init

   after the initial clone command without any directory path argument. To download
   only selected additional files, run one or more of the following commands.

The MIRTK group hosts Git submodule repositories with copies of the Boost_ and Eigen_
library header files required to build some of the core modules of the MIRTK
(see `MIRTK/Boost <https://github.com/MIRTK/Boost>`__ and
`MIRTK/Eigen <https://github.com/MIRTK/Eigen>`__ projects on GitHub).
Either install a compatible version of these libraries using the official installation
packages available for your operating system or initialize the desired Git submodules
before configuring the MIRTK build:

.. code-block:: bash

    git submodule update --init -- ThirdParty/Boost
    git submodule update --init -- ThirdParty/Eigen

The source files of core MIRTK modules are included in the top-level MIRTK repository
under the Modules/ subdirectory. Further optional packages which are developed and
distributed in their own respective Git repository are included as Git submodules
in the Packages/ directory. To also build these optional packages, initialize and
update the desired Git submodules before (re-)configuring the build files, e.g.:

.. code-block:: bash

    git submodule update --init -- Packages/Deformable
    git submodule update --init -- Packages/VolumetricMapping


.. _Git:   https://git-scm.com
.. _Boost: http://www.boost.org
.. _Eigen: http://eigen.tuxfamily.org


System Requirements
===================

**Operating System:**  Linux, Mac OS X, Microsoft Windows


Software License
================

The MIRTK is distributed under the terms of the
`Apache License Version 2 <http://www.apache.org/licenses/LICENSE-2.0>`__.
The license enables usage of MIRTK in both commercial and non-commercial applications,
without restrictions on the licensing applied to the combined work.

The MIRTK Git repository includes source files and references to Git submodule repositories
whose source files are covered by their own respective license terms which are compatible
with the MIRTK license. See the following links for details:

- `ThirdParty/Boost <https://github.com/MIRTK/Boost>`__: `Boost Software License Version 1.0 <http://www.boost.org/users/license.html>`__
- `ThirdParty/Eigen <https://github.com/MIRTK/Eigen>`__: `Mozilla Public License Version 2.0 <https://www.mozilla.org/en-US/MPL/2.0/>`__
- `ThirdParty/LBFGS <https://github.com/BioMedIA/MIRTK/tree/master/ThirdParty/LBFGS>`__: `The MIT License <https://opensource.org/licenses/MIT>`__
- `ThirdParty/NIfTI <https://github.com/BioMedIA/MIRTK/tree/master/ThirdParty/NIfTI>`__: `Public domain <https://en.wikipedia.org/wiki/Public_domain>`__
