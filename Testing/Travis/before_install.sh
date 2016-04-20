## Travis CI before_install script for Ubuntu 14.04 (Trusty) environment
set -e

# install pre-requisites on Ubuntu
if [ $TRAVIS_OS_NAME = linux ]; then

  sudo apt-get update -qq
  sudo apt-get install -y --no-install-recommends \
    freeglut3-dev \
    libarpack2-dev \
    libflann-dev \
    libgtest-dev \
    libnifti-dev \
    libpng12-dev \
    libsuitesparse-dev \
    libtbb-dev \
    libvtk6-dev

  mkdir $HOME/gtest-build && cd $HOME/gtest-build
  cmake /usr/src/gtest && make
  sudo mv -f libgtest.a libgtest_main.a /usr/lib

# install pre-requisites on OS X
elif [ $TRAVIS_OS_NAME = osx ]; then

  brew update
  brew tap homebrew/science
  brew install \
    arpack \
    flann \
    suite-sparse \
    tbb \
    vtk

  git clone --depth=1 https://github.com/google/googletest.git $HOME/gtest-source
  mkdir $HOME/gtest-build && cd $HOME/gtest-build
  cmake -DCMAKE_CXX_FLAGS=-std=c++11 -DBUILD_GMOCK=OFF -DBUILD_GTEST=ON ../gtest-source && make
  sudo make install

fi
