.. meta::
   :description: Get started with the Medical Image Registration ToolKit (MIRTK)
   :keywords:    image processing, image registration, IRTK, MIRTK, intro, quick start

===========
Get Started
===========

Install the software
--------------------

**Custom build**

Before you can start using the MIRTK or develop your own MIRTK applications,
you have to install the desired modules and available applications on your
system including their prerequisites. For a manual installation of the MIRTK,
please follow the :doc:`download` and :doc:`install` instructions which
list the third-party libraries required by each module and describe how to
build the MIRTK from its publicly available source code.


**Using Docker**

Instead of manually installing the MIRTK locally on your system, you can use
the pre-made `biomedia/mirtk <https://hub.docker.com/r/biomedia/mirtk/>`_
Docker image to run the :doc:`commands` inside a `Docker container`_.
For a guide to install and use Docker_, see the `official docs <https://docs.docker.com>`__.


**AppImage for Linux**

For Linux users, the easiest way to get started with MIRTK without the need for the Docker
runtime environment is the `MIRTK AppImage`_ available on Bintray_. An AppImage_ contains
all the required shared libraries and pre-built MIRTK commands. It can be executed on any
Linux system with a compatible minimum glibc version (>=2.15). With this AppImage, there is no
actual need for an installation. Simply download the file, make it executable, and copy it to
a directory that is in your PATH, e.g.,::

  wget -O mirtk https://bintray.com/schuhschuh/AppImages/download_file?file_path=MIRTK-latest-x86_64-glibc2.14.AppImage
  chmod a+x mirtk
  sudo mv mirtk /usr/bin

This AppImage is updated automatically when a change is committed to the master branch.


.. _BashCompletion:

Enable Bash completion
----------------------

**Custom build**

For information on how to enable auto-completion_ when running ``mirtk`` commands of
a manual MIRTK installation in the Bash_ shell, see the :doc:`install` instructions.


**Using Docker**

To enable Bash completion for running the MIRTK commands with Docker,
copy the `docker <https://raw.githubusercontent.com/docker/docker/master/contrib/completion/bash/docker>`__
and `docker-mirtk <https://raw.githubusercontent.com/BioMedIA/MIRTK/master/Docker/Completion/Bash/docker-mirtk>`__
completion scripts to ``/etc/bash_completion.d/`` on Linux or,
with the "bash-completion" Homebrew_ package installed on OS X,
to ``/usr/local/etc/bash_completion.d/``, respectively.

This can be done on Linux with the following Terminal commands::

  sudo curl -L https://raw.githubusercontent.com/docker/docker/master/contrib/completion/bash/docker       > /etc/bash_completion.d/docker
  sudo curl -L https://raw.githubusercontent.com/BioMedIA/MIRTK/master/Docker/Completion/Bash/docker-mirtk > /etc/bash_completion.d/docker-mirtk

On OS X with Homebrew, use these commands instead::

  brew install bash-completion
  curl -L https://raw.githubusercontent.com/docker/docker/master/contrib/completion/bash/docker        > $(brew --prefix)/etc/bash_completion.d/docker
  curl -L https://raw.githubusercontent.com/BioMedIA/MIRTK/master/Docker/Completion/Bash/docker-mirtk > $(brew --prefix)/etc/bash_completion.d/docker-mirtk

Alternatively, save the files to your home directory at, for example, ``$HOME/bash_completion/``
and add the following lines to your ``.bashrc`` (Linux) or ``.bash_profile`` (OS X) file::

  [ ! -f "$HOME/bash_completion/docker"       ] || . "$HOME/bash_completion/docker"
  [ ! -f "$HOME/bash_completion/docker-mirtk" ] || . "$HOME/bash_completion/docker-mirtk"


**Using AppImage**

To enable Bash completion for running the MIRTK commands in Bash using the AppImage for Linux,
copy `this file <https://raw.githubusercontent.com/BioMedIA/MIRTK/master/Scripts/mirtk_bash_completion.sh>`__
to ``/etc/bash_completion.d/``.

This can be done with the following Terminal command::

  sudo curl -L https://raw.githubusercontent.com/BioMedIA/MIRTK/master/Scripts/mirtk_bash_completion.sh > /etc/bash_completion.d/mirtk

Alternatively, save the file to your home directory at, for example, ``$HOME/bash_completion/``,
rename it to ``mirtk``, and add the following line to your ``.bashrc`` file::

  [ ! -f "$HOME/bash_completion/mirtk" ] || . "$HOME/bash_completion/mirtk"


Run the commands
----------------

The MIRTK installs a single executable named ``mirtk`` in the ``bin/`` directory
of the installation root directory. This executable is used to execute the various
MIRTK commands. The name of the command must be specified as first argument of the
``mirtk`` executable. Help about a command can be printed using the special ``help``
or ``help-rst`` commands followed by the name of the command for which help is requested::

    mirtk help <command>
    mirtk [-v] [-v] <command> [options]

When the ``-v`` option is given to ``mirtk`` before the command name, the path and
arguments of the command executable is printed before execution.

To run the commands inside a Docker_ container instead which does not require a local
installation of MIRTK, use the following command::

    docker run --rm --volume=<path>:/data biomedia/mirtk [-v] [-v] <command> [<options>]

This will download the `MIRTK Docker image`_ upon first execution. This image is
stored locally and will be reused for consecutive executions.
The ``--volume`` option of ``docker run`` mounts the specified directory path on the
host system to the ``/data`` directory inside the MIRTK Docker container.
The ``/data`` directory is the working directory of the MIRTK command.
The ``--rm`` option automatically deletes the MIRTK Docker container after the
command finished. Note that each ``docker run`` will create a new Docker container.
As these containers are meant to be used only once for each command execution, they
should be removed again after the command has finished.

For example, to print information about the NIfTI image file ``/path/to/my/images/image1.nii.gz``
using the :doc:`commands/info` command, execute the Docker command::

    docker run --rm --volume=/path/to/my/images:/data biomedia/mirtk info image1.nii.gz

See the :doc:`commands` page for a description of each command and the available options.

.. note::

   When Bash completion of MIRTK commands is enabled (see :ref:`BashCompletion`),
   press the <tab> key twice after typing "mirtk " or "docker run biomedia/mirtk "
   (incl. a space) to see a list of all commands. To only see a list of partial matches,
   type "mirtk eval" or  "docker run biomedia/mirtk eval", for example, and press <tab> twice.


Write your own application
--------------------------

For writing your own MIRTK command or an application which uses the MIRTK libraries,
we recommend a look at the source code of the applications included in the MIRTK.
The :doc:`API Reference <apidoc>` generated by Doxygen_ provides a more detailed
overview of the available MIRTK classes and their interfaces.

If you intend to contribute your applications in the future to the MIRTK distribution,
see the `code contribution <https://github.com/BioMedIA/MIRTK/blob/master/CONTRIBUTING.md>`__
guidelines for more information on how to contribute your code to the MIRTK source tree
or develop your own MIRTK Package.


.. _AppImage:           https://appimage.org/
.. _MIRTK AppImage:     https://bintray.com/schuhschuh/AppImages/MIRTK/master
.. _Bintray:            https://bintray.com/schuhschuh/AppImages/MIRTK/master
.. _Bash:               https://www.gnu.org/software/bash/
.. _auto-completion:    https://www.gnu.org/software/bash/manual/html_node/Programmable-Completion.html
.. _Homebrew:           http://brew.sh
.. _Doxygen:            http://www.doxygen.org/
.. _Docker:             http://www.docker.com
.. _Docker container:   https://www.docker.com/what-docker
.. _MIRTK Docker image: https://hub.docker.com/r/biomedia/mirtk/
