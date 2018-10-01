## Install dependencies on Ubuntu or OS X (using Homebrew)

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
WITH_FLTK=`norm_option_value "$WITH_FLTK" OFF`
WITH_CCACHE=`norm_option_value "$WITH_CCACHE" OFF`
BUILD_DEPS_WITH_CCACHE=`norm_option_value "$BUILD_DEPS_WITH_CCACHE" OFF`
FORCE_REBUILD_DEPS=`norm_option_value "$FORCE_REBUILD_DEPS" OFF`


# ------------------------------------------------------------------------------
# Auxiliary variables and functions
if [ $TRAVIS = ON ]; then
  os="$TRAVIS_OS_NAME"
else
  os=${1:-`uname`}
  travis_wait()
  {
    "$@"
  }
fi

run()
{
  echo "> $@"
  "$@"
  [ $? -eq 0 ] || exit 1
}


# ------------------------------------------------------------------------------
# install pre-requisites on Ubuntu
if [ $os = linux ] || [ $os = Linux ]; then
  cpu_cores=$(grep -c ^processor /proc/cpuinfo)

  deps=( \
    freeglut3-dev \
    libboost-math-dev \
    libboost-random-dev \
    libeigen3-dev \
    libnifti-dev \
    libpng12-dev \
  )

  [ $TESTING = OFF ] || deps=(${deps[@]} libgtest-dev)
  [ $WITH_TBB = OFF ] || deps=(${deps[@]} libtbb-dev)
  [ $WITH_FLANN = OFF ] || deps=(${deps[@]} libflann-dev)
  [ $WITH_ARPACK = OFF ] || deps=(${deps[@]} libarpack2-dev)

  if [ $WITH_UMFPACK = ON ]; then
    # see https://bugs.launchpad.net/ubuntu/+source/suitesparse/+bug/1333214
    # sudo add-apt-repository -y ppa:bzindovic/suitesparse-bugfix-1319687 || exit 1
    deps=(${deps[@]} libsuitesparse-dev)
  fi

  if [ $WITH_VTK = ON ]; then
    if [ -n "$LINUX_VTK_VERSION" ]; then
      VTK_VERSION="$LINUX_VTK_VERSION"
    fi
    if [ -z "$VTK_VERSION" ] || [ $VTK_VERSION = '6.0.0' ]; then
      deps=(${deps[@]} libvtk6-dev)
      VTK_VERSION=''
    fi
    if [ $WITH_FLTK = ON ]; then
      deps=(${deps[@]} libxi-dev libxmu-dev libxinerama-dev libxcursor-dev libcairo-dev libfltk1.3-dev)
    fi
  fi

  sudo apt-get update -qq || exit 1
  sudo apt-get install -y --no-install-recommends ${deps[@]} || exit 1

  if [ $TESTING = ON ]; then
    # libgtest-dev only install source files
    mkdir /tmp/gtest-build && cd /tmp/gtest-build
    run cmake /usr/src/gtest
    run make -j $cpu_cores
    run sudo mv -f libgtest.a libgtest_main.a /usr/lib
    cd || exit 1
    [ $TRAVIS = ON ] || rm -rf /tmp/gtest-build
  fi
fi


# ------------------------------------------------------------------------------
# install pre-requisites on OS X
if [ $os = osx ] || [ $os = Darwin ]; then
  cpu_cores=$(sysctl -n hw.ncpu)

  brew_install()
  {
    for dep in $@; do
      if $(brew ls $dep &> /dev/null); then
        brew unlink $dep && brew link $dep
        [ $? -eq 0 ] || exit 1
      else
        brew install $dep
      fi
    done
  }

  brew update > /dev/null || exit 1
  if [ $WITH_CCACHE = ON ]; then
    brew_install ccache
  fi
  brew_install eigen flann tbb
  if [ $WITH_ARPACK = ON ]; then
    brew_install arpack
  fi
  if [ $WITH_UMFPACK = ON ]; then
    brew_install suite-sparse
  fi
  if [ $WITH_VTK = ON ]; then
    VTK_VERSION="$MACOS_VTK_VERSION"
    if [ -z "$VTK_VERSION" ]; then
      echo "Installing VTK using Homebrew"
      brew_install vtk --without-python
    fi
    if [ $WITH_FLTK = ON ]; then
      brew_install fltk
    fi
  fi

  # download, build, and install gtest
  if [ $TESTING = ON ]; then
    run git clone --depth=1 https://github.com/google/googletest.git /tmp/gtest-source
    mkdir /tmp/gtest-build && cd /tmp/gtest-build
    [ $? -eq 0 ] || exit 1
    run cmake -DCMAKE_CXX_FLAGS=-std=c++11 -DBUILD_GMOCK=OFF -DBUILD_GTEST=ON ../gtest-source
    run make -j $cpu_cores
    run sudo make install
    cd || exit 1
    [ $TRAVIS = ON ] || rm -rf /tmp/gtest-build /tmp/gtest-source
  fi
fi


# ------------------------------------------------------------------------------
# Install specific VTK version from source
if [ $WITH_VTK = ON ] && [ -n "$VTK_VERSION" ]; then
  vtk_prefix="${VTK_PREFIX:-$HOME/VTK-$VTK_VERSION}"
  # build configuration
  cmake_args=(
    -DCMAKE_INSTALL_PREFIX="$vtk_prefix"
    -DCMAKE_BUILD_TYPE=Release
  )
  # pre-requisites to use system installations
  if [ $os = osx ] || [ $os = Darwin ]; then
    brew_install hdf5 netcdf jpeg libpng libtiff lz4
    cmake_args=("${cmake_args[@]}"
      -DVTK_USE_SYSTEM_HDF5=ON
      -DVTK_USE_SYSTEM_EXPAT=ON
      -DVTK_USE_SYSTEM_LIBXML2=ON
      -DVTK_USE_SYSTEM_ZLIB=ON
      -DVTK_USE_SYSTEM_NETCDF=ON
      -DVTK_USE_SYSTEM_JPEG=ON
      -DVTK_USE_SYSTEM_PNG=ON
      -DVTK_USE_SYSTEM_TIFF=ON
      -DVTK_USE_SYSTEM_LIBRARIES=ON
    )
  fi
  if [ $FORCE_REBUILD_DEPS = OFF ] && [ -d "$vtk_prefix/lib/cmake/vtk-${VTK_VERSION%.*}" ]; then
    # use previously cached VTK installation
    echo "Using cached VTK $VTK_VERSION installation in $vtk_prefix"
  else
    # whether to remove downloaded sources and build directory after installation
    if [ $TRAVIS = ON ]; then
      cleanup_vtk_build=OFF
    else
      cleanup_vtk_build=ON
    fi
    # custom build instead of Homebrew to take advantage of caching of minimal build
    cd /tmp
    echo "Downloading VTK $VTK_VERSION..."
    run curl -O "https://www.vtk.org/files/release/${VTK_VERSION%.*}/VTK-${VTK_VERSION}.tar.gz"
    run tar -xzf "VTK-${VTK_VERSION}.tar.gz"
    mkdir "VTK-${VTK_VERSION}/Build" && cd "VTK-${VTK_VERSION}/Build"
    [ $? -eq 0 ] || exit 1
    echo "Configuring VTK $VTK_VERSION..."
    cmake_args=("${cmake_args[@]}"
      -DVTK_Group_StandAlone=OFF
      -DVTK_Group_Rendering=OFF
      -DModule_vtkCommonCore=ON
      -DModule_vtkCommonDataModel=ON
      -DModule_vtkCommonExecutionModel=ON
      -DModule_vtkFiltersCore=ON
      -DModule_vtkFiltersHybrid=ON
      -DModule_vtkFiltersFlowPaths=ON
      -DModule_vtkFiltersGeneral=ON
      -DModule_vtkFiltersGeometry=ON
      -DModule_vtkFiltersParallel=ON
      -DModule_vtkFiltersModeling=ON
      -DModule_vtkImagingStencil=ON
      -DModule_vtkIOLegacy=ON
      -DModule_vtkIOXML=ON
      -DModule_vtkIOGeometry=ON
      -DModule_vtkIOPLY=ON
      -DModule_vtkIOXML=ON
    )
    if [ $WITH_CCACHE = ON ] && [ $BUILD_DEPS_WITH_CCACHE = ON ]; then
      cc_compiler=''
      cxx_compiler=''
      launcher=`which ccache`
      [ -z "$CC"  ] || cc_compiler=`which $CC`
      [ -z "$CXX" ] || cc_compiler=`which $CXX`
      if [ "$cc_compiler" = "${cc_compiler/ccache/}" ]; then
        echo "Using $launcher as C compiler launcher"
        cmake_args=("${cmake_args[@]}"
          -DCMAKE_C_COMPILER_LAUNCHER="$launcher"
        )
      fi
      if [ "$cxx_compiler" = "${cxx_compiler/ccache/}" ]; then
        echo "Using $launcher as C++ compiler launcher"
        cmake_args=("${cmake_args[@]}"
          -DCMAKE_CXX_COMPILER_LAUNCHER="$launcher"
        )
      fi
    fi
    run cmake "${cmake_args[@]}" ..
    echo "Configuring VTK $VTK_VERSION... done"
    echo "Building VTK $VTK_VERSION..."
    run make -j $cpu_cores
    echo "Building VTK $VTK_VERSION... done"
    echo "Installing VTK $VTK_VERSION..."
    run sudo make install
    echo "Installing VTK $VTK_VERSION... done"
    cd
    [ $TRAVIS = ON ] || rm -rf "/tmp/VTK-${VTK_VERSION}"
  fi
  mkdir -p "$HOME/.cmake/packages/VTK" || exit 1
  echo "$vtk_prefix/lib/cmake/vtk-${VTK_VERSION%.*}" > "$HOME/.cmake/packages/VTK/VTK-$VTK_VERSION"
fi
