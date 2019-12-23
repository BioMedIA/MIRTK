## Build Ubuntu Docker image for build of MIRTK
##
## The build of this image is time consuming. Moreover, this base image
## rarely changes and can thus be built and pushed to Docker Hub manually.
##
## How to build the image:
## - Change to the Docker/Ubuntu/ directory of the MIRTK source tree
## - Run "docker build -t biomedia/ubuntu:mirtk ."
##
## Upload image to Docker Hub:
## - Log in with "docker login" command if necessary
## - Push image using "docker push biomedia/ubuntu:mirtk"
##
## https://hub.docker.com/r/biomedia/ubuntu/

FROM ubuntu:18.04

LABEL Maintainer="Andreas Schuh <andreas.schuh.84@gmail.com>" \
      Description="Ubuntu with prerequisits for MIRTK installed" \
      Vendor="BioMedIA"

# No. of threads to use for build (--build-arg THREADS=8)
# By default, all available CPUs are used. When a Docker Machine is used,
# set the number of CPUs in the VirtualBox VM Settings.
ARG THREADS

# When VTK_VERSION='', the official libvtk7-dev package is used.
# Note, however, that this results in a Docker image that is about twice
# the size of the image when a custom VTK build without Qt, wrappers,
# and unused VTK modules is used instead!
ARG VTK_VERSION=8.2.0

ARG EIGEN_VERSION=3.3.7

ARG CXX_STANDARD=c++14

# Install prerequisites for MIRTK build
RUN NUM_CPUS=${THREADS:-`cat /proc/cpuinfo | grep processor | wc -l`} \
    && echo "Maximum number of build threads = $NUM_CPUS" \
    && apt-get update && apt-get install -y --no-install-recommends \
      software-properties-common \
    && apt-get update && apt-get install -y --no-install-recommends \
      wget \
      gcc \
      g++ \
      make \
      cmake \
      python3 \
      freeglut3-dev \
      libarpack2-dev \
      libboost-math-dev \
      libboost-random-dev \
      libflann-dev \
      libfftw3-dev \
      libgtest-dev \
      libinsighttoolkit4-dev \
      libnifti-dev \
      libpng-dev \
      libsuitesparse-dev \
      libtbb-dev \
      uuid-dev \
      zlib1g-dev \
    && mkdir /usr/src/gtest/build \
    && cd /usr/src/gtest/build \
    && cmake .. \
    && make -j $NUM_CPUS \
    && mv -f libgtest.a libgtest_main.a /usr/lib \
    && cd /usr/src \
    && rm -rf /usr/src/gtest/build \
    && \
    if [ -z ${VTK_VERSION} ]; then \
      apt-get install -y libeigen3-dev; \
    else \
      EIGEN_SOURCE_DIR=/usr/src/eigen-${EIGEN_VERSION} \
      && mkdir ${EIGEN_SOURCE_DIR} /usr/include/eigen3 \
      && cd ${EIGEN_SOURCE_DIR} \
      && wget -O archive.tar.bz2 http://bitbucket.org/eigen/eigen/get/${EIGEN_VERSION}.tar.bz2 \
      && tar vxjf archive.tar.bz2 --strip 1 \
      && mv signature_of_eigen3_matrix_library Eigen /usr/include/eigen3/ \
      && mv cmake/FindEigen3.cmake /usr/share/cmake-3.10/Modules/ \
      && cd /usr/src \
      && rm -rf ${EIGEN_SOURCE_DIR}; \
    fi \
    && \
    if [ -z ${VTK_VERSION} ]; then \
      apt-get install -y libvtk7-dev; \
    else \
      VTK_RELEASE=`echo ${VTK_VERSION} | sed s/\.[0-9]*$//` \
      && cd /usr/src \
      && wget http://www.vtk.org/files/release/${VTK_RELEASE}/VTK-${VTK_VERSION}.tar.gz \
      && tar -xvzf VTK-${VTK_VERSION}.tar.gz \
      && rm -f VTK-${VTK_VERSION}.tar.gz \
      && mkdir VTK-${VTK_VERSION}/Build \
      && cd VTK-${VTK_VERSION}/Build \
      && cmake \
        -D CMAKE_INSTALL_PREFIX=/usr/local \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_CXX_STANDARD=${CXX_STANDARD} \
        -D VTK_USE_SYSTEM_PNG=ON \
        -D VTK_USE_SYSTEM_ZLIB=ON \
        -D BUILD_SHARED_LIBS=ON \
        -D BUILD_EXAMPLES=OFF \
        -D BUILD_TESTING=OFF \
        -D BUILD_DOCUMENTATION=OFF \
        .. \
      && make -j $NUM_CPUS install \
      && cd /usr/src \
      && rm -rf /usr/src/VTK-${VTK_VERSION} \
      && ldconfig; \
    fi \
    && rm -rf /var/lib/apt/lists/*
