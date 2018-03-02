# MIRTK helloworld

`flip.cc` is just `Applications/src/flip-image.cc` renamed, and is supposed to
be an example program that you might try to build against 
[MIRTK](https://github.com/BioMedIA/MIRTK).

# Setting up MIRTK as a build dependency

Set `MIRTK_ROOT` to point to your MIRTK install area, and make sure the bin
directory is on your `PATH` and the lib area is on your library path. For
example, on Linux you might append these lines to your `.bashrc`:

```
export MIRTK_ROOT=/opt/mirtk
export PATH="$MIRTK_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$MIRTK_ROOT/lib:$LD_LIBRARY_PATH"
```

Then build this program with:

```
$ cd MIRTK/Examples/helloworld
$ mkdir build
$ cd build
$ cmake -D CMAKE_MODULE_PATH:PATH=$MIRTK_ROOT/lib/cmake/mirtk ..
$ make
$ make install
```

And run with:

```
$ ./flip
```
