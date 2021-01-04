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

CXX_STANDARD=11
TRAVIS=`norm_option_value "$TRAVIS" OFF`
TESTING=`norm_option_value "$TESTING" OFF`
WITH_ARPACK=`norm_option_value "$WITH_ARPACK" OFF`
WITH_UMFPACK=`norm_option_value "$WITH_UMFPACK" OFF`
WITH_ITK=`norm_option_value "$WITH_ITK" OFF`
WITH_VTK=`norm_option_value "$WITH_VTK" OFF`
WITH_TBB=`norm_option_value "$WITH_TBB" ON`
WITH_FLANN=`norm_option_value "$WITH_FLANN" ON`
WITH_FLTK=`norm_option_value "$WITH_FLTK" OFF`
WITH_CCACHE=`norm_option_value "$WITH_CCACHE" OFF`
BUILD_DEPS_WITH_CCACHE=`norm_option_value "$BUILD_DEPS_WITH_CCACHE" OFF`
FORCE_REBUILD_DEPS=`norm_option_value "$FORCE_REBUILD_DEPS" OFF`
DEBUG_VTK_BUILD=`norm_option_value "$DEBUG_VTK_BUILD" OFF`


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

  if [ -f /etc/lsb-release ]; then
    source /etc/lsb-release
  else
    echo "Script works on Ubuntu only" 1>&2
    exit 1
  fi
  if [ "$DISTRIB_ID" != "Ubuntu" ]; then
    echo "Script requires Ubuntu 14.04, 16.04, or 18.04" 1>&2
    exit 1
  fi

  cmake_cmd="$(which cmake)"
  if [ -n "$cmake_cmd" ]; then
    cmake_version="$("$cmake_cmd" --version | grep 'cmake version' | cut -d' ' -f3)"
    echo "Found CMake version $cmake_version"
    cmake_version_major="${cmake_version/.*}"
    [ $? -eq 0 -a -n "$cmake_version_major" ] || cmake_version_major=0
  else
    cmake_version_major=0
  fi
  if [ ${cmake_version_major} -lt 3 ]; then
    cmake_version=3.12.4
    echo "Installing CMake version $cmake_version"
    wget --quiet https://cmake.org/files/v${cmake_version%.*}/cmake-${cmake_version}-Linux-x86_64.tar.gz -O /tmp/cmake.tar.gz
    mkdir /opt/cmake-${cmake_version}
    tar xf /tmp/cmake.tar.gz -C /opt/cmake-${cmake_version} --strip-components=1 -h
    rm -f /tmp/cmake.tar.gz
    cmake_cmd="/opt/cmake-${cmake_version}/bin/cmake"
  fi

  deps=( \
    freeglut3-dev \
    libboost-math-dev \
    libboost-random-dev \
    libeigen3-dev \
    libnifti-dev \
    libpng-dev \
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
    if [ "$DISTRIB_CODENAME" = "trusty" ]; then
      if [ -z "$VTK_VERSION" ] || [ $VTK_VERSION = '6.0.0' ]; then
        deps=(${deps[@]} libvtk6-dev)
        VTK_VERSION=''
      fi
    elif [ "$DISTRIB_CODENAME" = "xenial" ]; then
      if [ -z "$VTK_VERSION" ] || [ $VTK_VERSION = '6.2.0' ]; then
        deps=(${deps[@]} libvtk6-dev python-vtk6)
        VTK_VERSION=''
      fi
    elif [ "$DISTRIB_CODENAME" = "bionic" ] || [ "$DISTRIB_CODENAME" = "focal" ] || [ "$DISTRIB_CODENAME" = "groovy" ]; then
      if [ $VTK_VERSION = '6.3.0' ]; then
        deps=(${deps[@]} libvtk6-dev)
        VTK_VERSION=''
      elif [ -z "$VTK_VERSION" ] || [ $VTK_VERSION = '7.1.1' ]; then
        deps=(${deps[@]} libvtk7-dev)
        VTK_VERSION=''
      fi
    elif [ -z "$VTK_VERSION" ]; then
      deps=(${deps[@]} libvtk7-dev)
    fi
    if [ $WITH_FLTK = ON ]; then
      deps=(${deps[@]} libxi-dev libxmu-dev libxinerama-dev libxcursor-dev libcairo-dev libfltk1.3-dev)
    fi
  fi

  if [ $WITH_ITK = ON ]; then
    deps=(${deps[@]} libinsighttoolkit4-dev libfftw3-dev uuid-dev)
  fi

  sudo apt-get update -qq || exit 1
  sudo apt-get install -y --no-install-recommends ${deps[@]} || exit 1

  if [ $TESTING = ON ]; then
    # libgtest-dev only install source files
    mkdir /tmp/gtest-build && cd /tmp/gtest-build
    run "$cmake_cmd" /usr/src/gtest
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
      brew_install vtk
    elif [ "$VTK_VERSION" = "8.2" ] || [ "$VTK_VERSION" == "9.0" ]; then
      echo "Installing VTK $VTK_VERSION using Homebrew"
      brew_install vtk@$VTK_VERSION
      VTK_VERSION=""  # skip installation from source code below
    fi
    if [ $WITH_FLTK = ON ]; then
      brew_install fltk
    fi
  fi
  if [ $WITH_ITK = ON ]; then
    brew_install itk fftw libuuid
  fi

  # download, build, and install gtest
  if [ $TESTING = ON ]; then
    run git clone --depth=1 https://github.com/google/googletest.git /tmp/gtest-source
    mkdir /tmp/gtest-build && cd /tmp/gtest-build
    [ $? -eq 0 ] || exit 1
    run "$cmake_cmd" -DCMAKE_CXX_STANDARD=$CXX_STANDARD -DBUILD_GMOCK=OFF -DBUILD_GTEST=ON ../gtest-source
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
    if [ ${VTK_VERSION/.*/} -lt 9 ]; then
      cmake_args+=(
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
    else
      cmake_args+=(
        -DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON
        -DVTK_MODULE_USE_EXTERNAL_VTK_expat=ON
        -DVTK_MODULE_USE_EXTERNAL_VTK_libxml2=ON
        -DVTK_MODULE_USE_EXTERNAL_VTK_zlib=ON
        -DVTK_MODULE_USE_EXTERNAL_VTK_netcdf=ON
        -DVTK_MODULE_USE_EXTERNAL_VTK_jpeg=ON
        -DVTK_MODULE_USE_EXTERNAL_VTK_png=ON
        -DVTK_MODULE_USE_EXTERNAL_VTK_tiff=ON
      )
    fi
  fi
  if [ $FORCE_REBUILD_DEPS = OFF ] && [ -d "$vtk_prefix/lib/cmake/vtk-${VTK_VERSION%.*}" ]; then
    # use previously cached VTK installation
    echo "Using cached VTK $VTK_VERSION installation in $vtk_prefix"
  else
    run rm -rf "$vtk_prefix"
    # custom build instead of Homebrew to take advantage of caching of minimal build
    cd /tmp
    echo "Downloading VTK $VTK_VERSION..."
    run curl -L -o "VTK-${VTK_VERSION}.tar.gz" "https://github.com/Kitware/VTK/archive/v${VTK_VERSION}.tar.gz"
    run tar -xzf "VTK-${VTK_VERSION}.tar.gz"
    mkdir "VTK-${VTK_VERSION}/Build" && cd "VTK-${VTK_VERSION}/Build"
    [ $? -eq 0 ] || exit 1
    [ "$DEBUG_VTK_BUILD" != "ON" ] || set -x
    echo "Configuring VTK $VTK_VERSION..."
    cmake_args+=(
      -DBUILD_TESTING=OFF
      -DCMAKE_CXX_STANDARD=$CXX_STANDARD
    )
    if [ ${VTK_VERSION/.*/} -lt 9 ]; then
      cmake_args+=(
        -DVTK_Group_Rendering=OFF
        -DVTK_Group_StandAlone=OFF
        -DVTK_WRAP_PYTHON=OFF
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
    else
      cmake_args+=(
        -DVTK_GROUP_ENABLE_Imaging=YES
        -DVTK_GROUP_ENABLE_MPI=DONT_WANT
        -DVTK_GROUP_ENABLE_Qt=DONT_WANT
        -DVTK_GROUP_ENABLE_Rendering=DONT_WANT
        -DVTK_GROUP_ENABLE_StandAlone=DONT_WANT
        -DVTK_GROUP_ENABLE_Views=DONT_WANT
        -DVTK_GROUP_ENABLE_Web=DONT_WANT
        -DVTK_MODULE_ENABLE_VTK_CommonCore=YES
        -DVTK_MODULE_ENABLE_VTK_CommonDataModel=YES
        -DVTK_MODULE_ENABLE_VTK_CommonExecutionModel=YES
        -DVTK_MODULE_ENABLE_VTK_FiltersCore=YES
        -DVTK_MODULE_ENABLE_VTK_FiltersHybrid=YES
        -DVTK_MODULE_ENABLE_VTK_FiltersFlowPaths=YES
        -DVTK_MODULE_ENABLE_VTK_FiltersGeneral=YES
        -DVTK_MODULE_ENABLE_VTK_FiltersGeometry=YES
        -DVTK_MODULE_ENABLE_VTK_FiltersParallel=YES
        -DVTK_MODULE_ENABLE_VTK_FiltersModeling=YES
        -DVTK_MODULE_ENABLE_VTK_ImagingStencil=YES
        -DVTK_MODULE_ENABLE_VTK_IOLegacy=YES
        -DVTK_MODULE_ENABLE_VTK_IOXML=YES
        -DVTK_MODULE_ENABLE_VTK_IOGeometry=YES
        -DVTK_MODULE_ENABLE_VTK_IOPLY=YES
        -DVTK_MODULE_ENABLE_VTK_IOXML=YES
      )
    fi
    if [ $WITH_CCACHE = ON ] && [ $BUILD_DEPS_WITH_CCACHE = ON ]; then
      cc_compiler=''
      cxx_compiler=''
      launcher=`which ccache`
      [ -z "$CC"  ] || cc_compiler=`which $CC`
      [ -z "$CXX" ] || cc_compiler=`which $CXX`
      if [ "$cc_compiler" = "${cc_compiler/ccache/}" ]; then
        echo "Using $launcher as C compiler launcher"
        cmake_args+=(-DCMAKE_C_COMPILER_LAUNCHER="$launcher")
      fi
      if [ "$cxx_compiler" = "${cxx_compiler/ccache/}" ]; then
        echo "Using $launcher as C++ compiler launcher"
        cmake_args+=(-DCMAKE_CXX_COMPILER_LAUNCHER="$launcher")
      fi
    fi
    run "$cmake_cmd" "${cmake_args[@]}" ..
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
  set +x
fi
