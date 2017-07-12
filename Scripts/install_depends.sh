## Install dependencies on Ubuntu or OS X (using Homebrew)
set -e

norm_option_value()
{
  if [ "$1" = on ] || [ "$1" = ON ] || [ "$1" = yes ] || [ "$1" = YES ] || [ "$1" = y ] || [ "$1" = Y ] || [ "$1" = 1 ] || [ "$1" = true ] || [ "$1" = TRUE ]; then
    echo ON
  elif [ "$1" = off ] || [ "$1" = OFF ] || [ "$1" = no ] || [ "$1" = NO ] || [ "$1" = n ] || [ "$1" = N ] || [ "$1" = 0 ] || [ "$1" = false ] || [ "$1" = FALSE ]; then
    echo OFF
  elif [ -n "$2" ]; then
    echo "$2"
  else
    echo OFF
  fi
}

TRAVIS=`norm_option_value "$TRAVIS" OFF`
TESTING=`norm_option_value "$TESTING" OFF`
WITH_ARPACK=`norm_option_value "$WITH_ARPACK" OFF`
WITH_UMFPACK=`norm_option_value "$WITH_UMFPACK" OFF`
WITH_VTK=`norm_option_value "$WITH_VTK" OFF`
WITH_TBB=`norm_option_value "$WITH_TBB" ON`
WITH_FLANN=`norm_option_value "$WITH_FLANN" ON`
WITH_CCACHE=`norm_option_value "$WITH_CCACHE" OFF`

if [ $TRAVIS = ON ]; then
  os="$TRAVIS_OS_NAME"
else
  os=${1:-`uname`}
  travis_wait()
  {
    "$@"
  }
fi

# install pre-requisites on Ubuntu
if [ $os = linux ] || [ $os = Linux ]; then

  deps=( \
    freeglut3-dev \
    libboost-math-dev \
    libboost-random-dev \
    libeigen3-dev \
    libnifti-dev \
    libpng12-dev \
  )
  [ $TESTING      = OFF ] || deps=(${deps[@]} libgtest-dev)
  [ $WITH_TBB     = OFF ] || deps=(${deps[@]} libtbb-dev)
  [ $WITH_FLANN   = OFF ] || deps=(${deps[@]} libflann-dev)
  [ $WITH_VTK     = OFF ] || deps=(${deps[@]} libvtk6-dev)
  [ $WITH_ARPACK  = OFF ] || deps=(${deps[@]} libarpack2-dev)
  [ $WITH_UMFPACK = OFF ] || {
    # see https://bugs.launchpad.net/ubuntu/+source/suitesparse/+bug/1333214
    sudo add-apt-repository -y ppa:bzindovic/suitesparse-bugfix-1319687
    deps=(${deps[@]} libsuitesparse-dev)
  }

  sudo apt-get update -qq
  sudo apt-get install -y --no-install-recommends ${deps[@]}

  if [ $TESTING = ON ]; then
    mkdir "$HOME/gtest-build" && cd "$HOME/gtest-build"
    cmake /usr/src/gtest && make
    sudo mv -f libgtest.a libgtest_main.a /usr/lib
    rm -rf "$HOME/gtest-build"
  fi

# install pre-requisites on OS X
elif [ $os = osx ] || [ $os = Darwin ]; then

  brew_install()
  {
    for dep in $@; do
      if $(brew ls --version $dep &> /dev/null) ; then
        brew unlink $dep && brew link $dep
      elif [[ $dep == arpack ]]; then
        travis_wait brew install $dep
      else
        brew install $dep
      fi
    done
  }

  brew update
  brew tap homebrew/science
  brew_install eigen flann tbb
  if [ $WITH_CCACHE = ON ]; then
    brew_install ccache
  fi
  if [ $WITH_ARPACK = ON ]; then
    brew_install arpack
  fi
  if [ $WITH_UMFPACK = ON ]; then
    brew_install suite-sparse
  fi
  if [ $WITH_VTK = ON ]; then
    brew_install vtk
  fi
  if [ $TESTING = ON ]; then
    cleanup_gtest_source='OFF'
    if [ -d "$HOME/gtest-source" ]; then
      if [ $TRAVIS = ON ]; then
        # directory has been restored from cache
        cd "$HOME/gtest-source" && git update
      fi
    else
      git clone --depth=1 https://github.com/google/googletest.git "$HOME/gtest-source"
      cleanup_gtest_source='ON'
    fi
    mkdir "$HOME/gtest-build" && cd "$HOME/gtest-build"
    cmake -DCMAKE_CXX_FLAGS=-std=c++11 -DBUILD_GMOCK=OFF -DBUILD_GTEST=ON ../gtest-source && make
    sudo make install
    if [ $cleanup_gtest_source = ON ]; then
      rm -rf "$HOME/gtest-source"
    fi
    rm -rf "$HOME/gtest-build"
  fi

fi
