## Build of Docker image for execution of MIRTK commands within a Docker
## container with all modules and their prerequisites available in the image

FROM ubuntu:14.04

MAINTAINER Andreas Schuh <andreas.schuh.84@gmail.com>
LABEL Description="Medical Image Registration ToolKit (MIRTK)" Vendor="BioMedIA"

RUN apt-get update && apt-get install -y --no-install-recommends \
      wget \
      gcc \
      g++ \
      make \
      cmake \
      python \
      freeglut3-dev \
      libarpack2-dev \
      libflann-dev \
      libgtest-dev \
      libnifti-dev \
      libpng-dev \
      libsuitesparse-dev \
      libtbb-dev \
      zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /usr/src/gtest/build \
    && cd /usr/src/gtest/build \
    && cmake .. \
    && make \
    && mv -f libgtest.a libgtest_main.a /usr/lib \
    && cd /usr/src \
    && rm -rf /usr/src/gtest/build

RUN cd /usr/src \
    && wget http://www.vtk.org/files/release/6.3/VTK-6.3.0.tar.gz \
    && tar -xzf VTK-6.3.0.tar.gz \
    && rm -f VTK-6.3.0.tar.gz \
    && mkdir VTK-6.3.0/Build \
    && cd VTK-6.3.0/Build \
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
    && rm -rf /usr/src/VTK-6.3.0 \
    && ldconfig

COPY . /usr/src/MIRTK
RUN mkdir /usr/src/MIRTK/Build \
    && cd /usr/src/MIRTK/Build \
    && cmake \
      -D CMAKE_INSTALL_PREFIX=/usr/local \
      -D CMAKE_BUILD_TYPE=Release \
      -D BUILD_ALL_MODULES=ON \
      -D BUILD_SHARED_LIBS=ON \
      -D BUILD_APPLICATIONS=ON \
      -D BUILD_TESTING=OFF \
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
    && make install \
    && cd /usr/src \
    && rm -rf /usr/src/MIRTK

ENTRYPOINT ["/usr/bin/python", "/usr/local/bin/mirtk"]
CMD ["help"]

WORKDIR /data
