## Build Docker image for execution of MIRTK commands within a Docker
## container with all modules and applications available in the image
##
## How to build the image:
## - Change to top-level directory of MIRTK source tree
## - Run "docker build --build-arg VCS_REF=`git rev-parse --short HEAD` -t <user>/mirtk:latest ."
##
## Upload image to Docker Hub:
## - Log in with "docker login" if necessary
## - Push image using "docker push <user>/mirtk:latest"
##
## Note: The "biomedia/mirtk" repository on Docker Hub is automatically
##       updated with a new Docker image build when a change is committed
##       to the (develop/master) branch of the MIRTK GitHub repository.
##       Pushing manually build images to biomedia/mirtk is not possible
##       (unless the automatic build has been replaced by a manual build).
##
## https://hub.docker.com/r/biomedia/mirtk/

# Pre-made Ubuntu system with all MIRTK dependencies installed
FROM biomedia/ubuntu:18.04-mirtk

LABEL Maintainer="Andreas Schuh <andreas.schuh.84@gmail.com>" \
      Description="Medical Image Registration ToolKit (MIRTK)" \
      Vendor="BioMedIA"

# Git repository and commit SHA from which this Docker image was built
# (see https://microbadger.com/#/labels)
ARG VCS_REF
LABEL org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/BioMedIA/MIRTK"

# No. of threads to use for build (--build-arg THREADS=8)
# By default, all available CPUs are used. When a Docker Machine is used,
# set the number of CPUs in the VirtualBox VM Settings.
ARG THREADS

# Whether to build and run MIRTK tests before installation
# Override default with docker build --build-arg BUILD_TESTING=ON
ARG BUILD_TESTING=OFF

# Build and install MIRTK
COPY . /usr/src/mirtk
RUN ls /usr/src/mirtk \
    && NUM_CPUS=${THREADS:-`cat /proc/cpuinfo | grep processor | wc -l`} \
    && echo "Maximum number of build threads = $NUM_CPUS" \
    && mkdir /usr/src/mirtk/Build \
    && cd /usr/src/mirtk/Build \
    && cmake \
      -D CMAKE_INSTALL_PREFIX=/usr/local \
      -D CMAKE_BUILD_TYPE=Release \
      -D BUILD_SHARED_LIBS=ON \
      -D BUILD_APPLICATIONS=ON \
      -D BUILD_TESTING=${BUILD_TESTING} \
      -D BUILD_DOCUMENTATION=OFF \
      -D BUILD_CHANGELOG=OFF \
      -D MODULE_Common=ON \
      -D MODULE_Numerics=ON \
      -D MODULE_Image=ON \
      -D MODULE_IO=ON \
      -D MODULE_PointSet=ON \
      -D MODULE_Transformation=ON \
      -D MODULE_Registration=ON \
      -D MODULE_Deformable=ON \
      -D MODULE_DrawEM=ON \
      -D MODULE_Mapping=ON \
      -D MODULE_Scripting=ON \
      -D MODULE_Viewer=OFF \
      -D WITH_ARPACK=ON \
      -D WITH_FLANN=ON \
      -D WITH_MATLAB=OFF \
      -D WITH_NiftiCLib=ON \
      -D WITH_PNG=ON \
      -D WITH_PROFILING=ON \
      -D WITH_TBB=ON \
      -D WITH_UMFPACK=ON \
      -D WITH_ITK=ON \
      -D WITH_VTK=ON \
      -D WITH_ZLIB=ON \
      .. \
    && if [ ${BUILD_TESTING} = ON ]; then make -j $NUM_CPUS && make test; fi \
    && make -j $NUM_CPUS install \
    && cd /usr/src \
    && rm -rf /usr/src/mirtk \
    && echo "/usr/local/lib/mirtk" > /etc/ld.so.conf.d/mirtk.conf \
    && ldconfig

# Make "mirtk" the default executable for application containers
ENTRYPOINT ["/usr/bin/python3", "/usr/local/bin/mirtk"]
CMD ["help"]

# Assume user data volume to be mounted at /data
#   docker run --volume=/path/to/data:/data
WORKDIR /data
