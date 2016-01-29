===================
CMake BASIS Modules
===================

This directory contains the CMake modules of the [CMake BASIS][1] project only.
These modules are required by any project which takes advantage of the extended
CMake commands of CMake BASIS. Other components of CMake BASIS such as the
CMake BASIS Utilities (a library of common functions for each supported programming
language) and CMake BASIS Tools (e.g., the ```basisproject``` tool) are installed
through the CMake BASIS project.

License
=======

Copyright (c) 2011-2013 University of Pennsylvania   <br />
Copyright (c) 2013-2015 Andreas Schuh                <br />
Copyright (c) 2013-2014 Carnegie Mellon University

CMake BASIS is available under a BSD compatible license. The complete license text
can be found in the [COPYING.txt](/COPYING.txt) file.

Installation
============

Developers requiring only the CMake BASIS Modules are encouraged to include the
[CMake BASIS Modules][2] files directly in their own Git controlled project source
tree, either as Git subtree or (shallow) submodule.

To utilize these modules, we recommend the use of the ```basis-modules``` project
template instead of the ```basis``` template which requires a separate build of
the entire [CMake BASIS][3] project.

**TODO**: Add template files to CMake BASIS and link them here.


Using CMake modules as subtree
------------------------------

To add the CMake BASIS Modules as subtree to your project under the subdirectory path ```basis/```,
use the following two commands. The first adds a new remote which simplifies the following
commands:

```bash
git remote add -f basis-modules https://github.com/schuhschuh/cmake-basis-modules.git
git subtree add --prefix basis/ basis-modules master --squash
```

In order to update the modules at a later date to incorporate changes of the CMake
BASIS Modules into your project, use the following commands:


```bash
git fetch basis-modules master
git subtree pull --prefix basis/ basis-modules master --squash
```


Adding CMake modules as submodule
---------------------------------

An alternative to the ```git subtree``` command to add the CMake BASIS Modules to
your project, you can use ```git submodule``` instead. For a comparison of the two
commands and their ups and downs, read some of the many tutorials available online.

To add the CMake BASIS Modules as submodule to your project under the subdirectory path ```basis/```,
use the following commands. The ```.gitmodules``` file which records the added submodules
and the URL of the remote repository must be committed to your project.

```bash
git submodule add --depth=1 https://github.com/schuhschuh/cmake-basis-modules.git basis
git add .gitmodules
git commit -m 'add: CMake BASIS Modules'
```

In order to update the modules at a later date to incorporate changes of the CMake
BASIS Modules into your project, use the following commands:

```bash
cd basis/                                        # change to submodule directory
git checkout master                              # checkout master branch
git pull                                         # pull remote changes
cd ..                                            # change back to main repository
git add basis                                    # stash change of submodule SHA
git commit -m 'mod: Update CMake BASIS Modules'  # commit submodule change
```

First, you have to pull the changes from the remote repository and merge them into
your local submodule repository. Then you have to update the submodule commit SHA
recorded in your main repository to point to the latest commit of the CMake BASIS
Modules.


[1]: http://opensource.andreasschuh.com/cmake-basis
[2]: https://github.com/schuhschuh/cmake-basis-modules
[3]: https://github.com/schuhschuh/cmake-basis
