# Contributing to MIRTK

This file documents ways to contribute to the MIRTK.

## Reporting issues

The MIRTK uses Github's [issue tracker](https://github.com/BioMedIA/MIRTK/issues?state=open)
to record and track the list of past and present issues. Before submitting a new issue,
please consider browsing this list to check whether your problem is being or has already
been dealt with. Together with a meaningful title, please provide detailed instructions to
help reproduce the issue. Important pieces of information include the OS, build context,
execution steps and expected outcome. 

Technical support may also be requested by speaking with the team directly via
[Gitter](https://gitter.im/BioMedIA/MIRTK). 

## Code contributions

The MIRTK maintainers welcome code contributions from MIRTK users, preferrably via
[pull requests](https://github.com/BioMedIA/MIRTK/pulls) against the master branch.
Please provide as much information as possible to help with the code review, such as why
this change is needed, what it does and how to check whether the code works as intended.

If you intend to contribute your applications to the MIRTK distribution, develop your tools
either in a fork of the MIRTK Git repository or your own MIRTK Package in a separate Git repository
with the common file organization and the appropriate [CMake](https://cmake.org) configuration
files. See some of the already available modules in the `Packages/` directory of the
MIRTK source tree as reference. Examples of MIRTK Packages are the
[Deformable](https://github.com/MIRTK/Deformable) and
[Volumetric Mapping](https://github.com/MIRTK/VolumetricMapping) modules,
each developed in their own respective GitHub repository.
Template files for the creation of a new MIRTK Package can be found in the
[Templates](Templates/mirtk-module) directory. These template files can be customized using
the [basisproject](https://cmake-basis.github.io/howto/use-and-customize-templates.html#use-a-template)
tool of [CMake BASIS](https://cmake-basis.github.io), but can also
simply be copied and edited by hand. But do not copy the `_conf.py` file from the
`Templates/mirtk-module/*/` directory.
