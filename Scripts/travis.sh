## Travis CI script
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

run()
{
  echo "> $@"
  $@
}

# option CIBUILD to allow disabling of build step
# in order for installed dependencies to be cached
CIBUILD=`norm_option_value "$CIBUILD" ON`
if [ $CIBUILD = OFF ]; then
  exit 0
fi

TRAVIS=`norm_option_value "$TRAVIS" OFF`
TESTING=`norm_option_value "$TESTING" OFF`
WITH_ARPACK=`norm_option_value "$WITH_ARPACK" OFF`
WITH_UMFPACK=`norm_option_value "$WITH_UMFPACK" OFF`
WITH_VTK=`norm_option_value "$WITH_VTK" OFF`
WITH_TBB=`norm_option_value "$WITH_TBB" ON`
WITH_FLANN=`norm_option_value "$WITH_FLANN" ON`
WITH_FLTK=`norm_option_value "$WITH_FLTK" OFF`
WITH_CCACHE=`norm_option_value "$WITH_CCACHE" OFF`

if [ $TRAVIS = ON ]; then
  os="$TRAVIS_OS_NAME"
else
  os=${1:-`uname`}
fi
if [ $os = linux ] || [ $os = Linux ]; then
  cpu_cores=$(grep -c ^processor /proc/cpuinfo)
elif [ $os = osx ] || [ $os = Darwin ]; then
  cpu_cores=$(sysctl -n hw.ncpu)
else
  cpu_cores=1
fi

modules=(Common Numerics Image IO Transformation Registration DrawEM)
if [ $WITH_VTK = ON ]; then
  modules=(${modules[@]} PointSet Deformable Mapping)
  if [ $WITH_FLTK = ON ]; then
    modules=(${modules[@]} Viewer)
  fi
fi

cmake_args=
for module in ${modules[@]}; do
  cmake_args=(${cmake_args[@]} -D MODULE_${module}=ON)
done

mkdir Build && cd Build
run cmake \
      -D CMAKE_INSTALL_PREFIX=/usr \
      -D CMAKE_BUILD_TYPE=Release \
      -D BUILD_SHARED_LIBS=ON \
      -D BUILD_APPLICATIONS=ON \
      -D BUILD_TESTING=$TESTING \
      -D BUILD_DOCUMENTATION=OFF \
      -D BUILD_CHANGELOG=OFF \
      -D WITH_CCACHE=$WITH_CCACHE \
      -D WITH_ARPACK=$WITH_ARPACK \
      -D WITH_FLANN=$WITH_FLANN \
      -D WITH_MATLAB=OFF \
      -D WITH_NiftiCLib=ON \
      -D WITH_PNG=ON \
      -D WITH_PROFILING=ON \
      -D WITH_TBB=$WITH_TBB \
      -D WITH_UMFPACK=$WITH_UMFPACK \
      -D WITH_VTK=$WITH_VTK \
      -D WITH_ZLIB=ON \
      ${cmake_args[@]} \
      ..

run make -j $cpu_cores
[ $TESTING = OFF ] || make test
