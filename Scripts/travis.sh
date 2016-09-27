## Travis CI script
set -e

norm_option_value()
{
  if [ $1 = on ] || [ $1 = ON ] || [ $1 = yes ] || [ $1 = YES ] || [ $1 = y ] || [ $1 = Y ] || [ $1 = 1 ]; then
    echo ON
  else
    echo OFF
  fi
}

TESTING=`norm_option_value "$TESTING"`
WITH_VTK=`norm_option_value "$WITH_VTK"`

modules=(Common Numerics Image IO Transformation Registration DrawEM)
if [ $WITH_VTK = ON ]; then
  modules=(${modules[@]} PointSet Deformable Mapping)
fi

cmake_args=
for module in ${modules[@]}; do
  cmake_args=(${config[@]} -D MODULE_${module}=ON)
done

mkdir Build && cd Build
cmake -D CMAKE_INSTALL_PREFIX=$HOME/local \
      -D CMAKE_BUILD_TYPE=Release \
      -D BUILD_SHARED_LIBS=ON \
      -D BUILD_APPLICATIONS=ON \
      -D BUILD_TESTING=$TESTING \
      -D BUILD_DOCUMENTATION=OFF \
      -D BUILD_CHANGELOG=OFF \
      -D WITH_ARPACK=ON \
      -D WITH_FLANN=ON \
      -D WITH_MATLAB=OFF \
      -D WITH_NiftiCLib=ON \
      -D WITH_PNG=ON \
      -D WITH_PROFILING=ON \
      -D WITH_TBB=ON \
      -D WITH_UMFPACK=ON \
      -D WITH_VTK=$WITH_VTK \
      -D WITH_ZLIB=ON \
      ${cmake_args[@]} \
      ..

make -j 8
[ $TESTING = OFF ] || make test
