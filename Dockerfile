## Build Docker image for execution of MIRTK commands within a Docker
## container with all modules and applications available in the image

# Pre-made Ubuntu system with all MIRTK dependencies installed
FROM biomedia/ubuntu:mirtk

MAINTAINER Andreas Schuh <andreas.schuh.84@gmail.com>
LABEL Description="Medical Image Registration ToolKit (MIRTK)" Vendor="BioMedIA"

# Whether to build and run MIRTK tests before installation
# Override default with docker build --build-arg BUILD_TESTING=ON
ARG BUILD_TESTING=OFF

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
