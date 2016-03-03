The MIRTK installs a single executable script file named ``mirtk`` in the
``bin/`` directory of the installation root directory. This executable is
used to execute the various MIRTK commands. The name of the command must
be specified as first argument of the ``mirtk`` executable. Help about
a command can be printed using the special ``help`` or ``help-rst`` commands
followed by the name of the command for which help is requested.

Usage::

    mirtk help <command>
    mirtk [-v] [-v] <command> [options]

When the ``-v`` option is given before the command name, the path and arguments
of the command executable is printed before execution.

.. note::

   The MIRTK installation includes a `Bash <https://www.gnu.org/software/bash/>`__
   completions script which enables auto-completion of the available MIRTK commands.
   See :ref:`InstallationSteps` on how to activate auto-completion. When enabled, press the
   <tab> key twice after typing "mirtk " (incl. a space) to see a list of all commands.
   To only see a list of partial matches type, for example, "mirtk eval" and press <tab> twice.
