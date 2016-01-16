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
[pull requests](https://github.com/BioMedIA/MIRTK/pulls) against the develop branch.
Please provide as much information as possible to help with the code review, such as why
this change is needed, what it does and how to check whether the code works as intended.

If you have developed your own MIRTK Package and would like it to be added to the MIRTK
distribution, please send a message to the mailing list or chat with the MIRTK team
on [Gitter](https://gitter.im/BioMedIA/MIRTK). Examples of MIRTK Packages
are the [Deformable](https://github.com/MIRTK/Deformable) and
[Volumetric Mapping](https://github.com/MIRTK/VolumetricMapping) modules,
each developed in their own respective GitHub repository. Template files for the creation
of a new MIRTK Package can be found in the [Templates](Templates/mirtk-module) directory.
These template files can be customized using the
[basisproject](http://opensource.andreasschuh.com/cmake-basis/howto/use-and-customize-templates.html#use-a-template)
tool of [CMake BASIS](http://opensource.andreasschuh.com/cmake-basis/),
but can also simply be copied and edited by hand.
