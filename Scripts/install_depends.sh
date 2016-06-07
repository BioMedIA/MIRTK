## Install dependencies on Ubuntu or OS X (using Homebrew)
set -e

if [ $TRAVIS = true ]; then
  os="$TRAVIS_OS_NAME"
else
  os=${1:-`uname`}
fi

# install pre-requisites on Ubuntu
if [ $os = linux ] || [ $os = Linux ]; then

  # see https://bugs.launchpad.net/ubuntu/+source/suitesparse/+bug/1333214
  sudo add-apt-repository -y ppa:bzindovic/suitesparse-bugfix-1319687

  sudo apt-get update -qq
  sudo apt-get install -y --no-install-recommends \
    freeglut3-dev \
    libarpack2-dev \
    libboost-math-dev \
    libboost-random-dev \
    libeigen3-dev \
    libflann-dev \
    libgtest-dev \
    libnifti-dev \
    libpng12-dev \
    libsuitesparse-dev \
    libtbb-dev \
    libvtk6-dev

  mkdir "$HOME/gtest-build" && cd "$HOME/gtest-build"
  cmake /usr/src/gtest && make
  sudo mv -f libgtest.a libgtest_main.a /usr/lib
  rm -rf "$HOME/gtest-build"

# install pre-requisites on OS X
elif [ $os = osx ] || [ $os = Darwin ]; then

  brew update
  brew tap homebrew/science
  brew install \
    arpack \
    eigen \
    flann \
    suite-sparse \
    tbb \
    vtk

  git clone --depth=1 https://github.com/google/googletest.git "$HOME/gtest-source"
  mkdir "$HOME/gtest-build" && cd "$HOME/gtest-build"
  cmake -DCMAKE_CXX_FLAGS=-std=c++11 -DBUILD_GMOCK=OFF -DBUILD_GTEST=ON ../gtest-source && make
  sudo make install
  rm -rf "$HOME/gtest-source" "$HOME/gtest-build"

fi
