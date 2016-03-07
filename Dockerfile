## Build of Docker image for execution of MIRTK commands within a Docker
## container with all MIRTK modules and applications available in the image
FROM ubuntu:14.04

MAINTAINER Andreas Schuh <andreas.schuh.84@gmail.com>
LABEL Description="Medical Image Registration ToolKit (MIRTK)" Vendor="BioMedIA"

# When no VTK_VERSION is set, the official libvtk6-dev package is used.
# Note, however, that this results in a Docker image that is about twice
# the size of the image when a custom VTK build without Qt, wrappers,
# and unused VTK modules is used instead! It is recommended to build
# the Docker image with "--build-arg VTK_VERSION=6.3.0" option.
ARG VTK_VERSION

# Whether to build and run MIRTK tests before installation ("--build-arg BUILD_TESTING=ON")
ARG BUILD_TESTING=OFF

# Install prerequisites
RUN apt-get update && apt-get install -y --no-install-recommends \
      gcc \
      g++ \
      make \
      cmake \
      python \
      freeglut3-dev \
      libarpack2-dev \
      libflann-dev \
      libnifti-dev \
      libpng-dev \
      libsuitesparse-dev \
      libtbb-dev \
      zlib1g-dev \
    && \
    if [ ${BUILD_TESTING} = ON ]; then \
      apt-get install -y libgtest-dev \
      && mkdir /usr/src/gtest/build \
      && cd /usr/src/gtest/build \
      && cmake .. \
      && make \
      && mv -f libgtest.a libgtest_main.a /usr/lib \
      && cd /usr/src \
      && rm -rf /usr/src/gtest/build; \
    fi \
    && \
    if [ -z ${VTK_VERSION} ]; then \
      apt-get install -y libvtk6-dev; \
    else \
      VTK_RELEASE=`echo ${VTK_VERSION} | sed s/\.[0-9]*$//` \
      && apt-get install -y wget \
      && cd /usr/src \
      && wget http://www.vtk.org/files/release/${VTK_RELEASE}/VTK-${VTK_VERSION}.tar.gz \
      && tar -xvzf VTK-${VTK_VERSION}.tar.gz \
      && rm -f VTK-${VTK_VERSION}.tar.gz \
      && mkdir VTK-${VTK_VERSION}/Build \
      && cd VTK-${VTK_VERSION}/Build \
      && cmake \
        -D CMAKE_INSTALL_PREFIX=/usr/local \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_CXX_FLAGS=-std=c++11 \
        -D VTK_USE_SYSTEM_PNG=ON \
        -D VTK_USE_SYSTEM_ZLIB=ON \
        -D BUILD_SHARED_LIBS=ON \
        -D BUILD_EXAMPLES=OFF \
        -D BUILD_TESTING=OFF \
        -D BUILD_DOCUMENTATION=OFF \
        .. \
      && make install \
      && cd /usr/src \
      && rm -rf /usr/src/VTK-${VTK_VERSION} \
      && ldconfig; \
    fi \
    && rm -rf /var/lib/apt/lists/*

# Build and install MIRTK
COPY . /usr/src/MIRTK
RUN mkdir /usr/src/MIRTK/Build \
    && cd /usr/src/MIRTK/Build \
    && cmake \
      -D CMAKE_INSTALL_PREFIX=/usr/local \
      -D CMAKE_BUILD_TYPE=Release \
      -D BUILD_ALL_MODULES=ON \
      -D BUILD_SHARED_LIBS=ON \
      -D BUILD_APPLICATIONS=ON \
      -D BUILD_TESTING=${BUILD_TESTING} \
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
      -D WITH_VTK=ON \
      -D WITH_ZLIB=ON \
      .. \
    && if [ ${BUILD_TESTING} = ON ]; then make && make test; fi \
    && make install \
    && cd /usr/src \
    && rm -rf /usr/src/MIRTK

# Make "mirtk" the default executable for application containers
ENTRYPOINT ["python", "/usr/local/bin/mirtk"]
CMD ["help"]

# Assume user data volume to be mounted at /data
#   docker run --volume=/path/to/data:/data
WORKDIR /data
