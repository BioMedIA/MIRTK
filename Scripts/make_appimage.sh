#!/bin/bash
set -e

########################################################################
# Package the binaries built on Travis CI as an AppImage
#
# For more information, see http://appimage.org/
########################################################################

# Output something before any command can fail (cf. 'set -e')
# If this is the only output in build log, use 'set -ex' above
echo "Creating AppImage"
echo

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
WITH_VTK=`norm_option_value "$WITH_VTK" OFF`
TRANSFER=`norm_option_value "$TRANSFER_AppImage" OFF`
AppImage_VIEWER=`norm_option_value "$AppImage_VIEWER" OFF`
AppImage_LATEST=`norm_option_value "$AppImage_LATEST" OFF`
PYTHON="/usr/bin/env python"  # target system Python interpreter

if [ $TRAVIS = ON ]; then
  sudo apt-get install -y tree
  cd "$TRAVIS_BUILD_DIR/Build"
  OUTDIR="$TRAVIS_BUILD_DIR/Deploy"
else
  OUTDIR="$(cd "$(dirname "$BASH_SOURCE")/.." && pwd)/Deploy"
fi

APP='MIRTK'
LOWERAPP=${APP,,}
RELEASE=''
if [ $TRAVIS = ON ]; then
  if [ -n "$TRAVIS_TAG" ]; then
    VERSION="$TRAVIS_TAG"
    RELEASE="$TRAVIS_TAG"
  else
    VERSION="$TRAVIS_COMMIT"
  fi
else
  VERSION=$(git describe --tags --exact-match --abbrev=0 HEAD 2> /dev/null || true)
  [ -n "$VERSION" ] || VERSION="$(git rev-parse --short HEAD)"
fi
if [ -z "$RELEASE" ]; then
  if [ $AppImage_LATEST = ON ]; then
    RELEASE="latest"  # overwrite previous "latest" BinTray AppImage
  else
    RELEASE="$VERSION"
  fi
fi
APPDIR="/tmp/$APP/$APP.AppDir"
export ARCH=$(arch)

echo "Name:    $APP"
echo "Version: $VERSION"
echo "Release: $RELEASE"
echo "Arch:    $ARCH"
echo
echo "Populating AppDir..."
make install DESTDIR="$APPDIR"

cd "$APPDIR/.."
wget -q https://github.com/probonopd/AppImages/raw/master/functions.sh -O ./functions.sh
. ./functions.sh
cd "$APPDIR"

########################################################################
# Include pre-built binary MIRTK view command (source not included)
########################################################################

if [ $AppImage_VIEWER = ON ]; then
  sudo apt-get install -y libgsl0ldbl
  wget -O view https://bintray.com/schuhschuh/generic/download_file?file_path=ubuntu-14.04%2Fview
  chmod a+x view
  if [ -d usr/lib/mirtk/tools ]; then
    mv view usr/lib/mirtk/tools/
  else
    mv view usr/lib/tools/
  fi
fi

########################################################################
# Copy in dependencies that may not be available on all target systems
########################################################################

if [ $WITH_VTK = ON ] && [ -n "$VTK_VERSION" ]; then
  export LD_LIBRARY_PATH="${VTK_PREFIX:-$HOME/VTK-$VTK_VERSION}/lib:$LD_LIBRARY_PATH"
fi
copy_deps
move_lib

########################################################################
# Delete stuff that should not go into the AppImage
########################################################################

# Delete dangerous libraries
# See https://github.com/probonopd/AppImages/blob/master/excludelist
delete_blacklisted

if [ -d usr/lib/$ARCH-linux-gnu ]; then
  mv usr/lib/$ARCH-linux-gnu/* usr/lib/
  rmdir usr/lib/$ARCH-linux-gnu
fi

if [ -d usr/lib/mirtk ]; then
  mv usr/lib/mirtk/* usr/lib/
  rmdir usr/lib/mirtk
  sed -i s:lib/mirtk:lib:g usr/bin/mirtk
fi

rm -f usr/bin/uninstall-mirtk || true
rm -rf usr/include || true
rm -rf usr/lib/cmake || true
rm -rf usr/doc || true
rm -rf usr/share || true
rmdir usr/lib64 || true
rmdir usr/lib/mesa || true

find . -name *.so -or -name *.so.* -exec strip {} \;
for f in $(find . -type f -executable -exec file -- {} \; | grep ELF | cut -d: -f1); do
  strip "$f"
done

########################################################################
# Fix she-bang interpreter specifications
########################################################################

for f in $(grep -rl '^#!.*python' usr); do
  echo -n "$f: replace '$(head -n1 "$f")'"
  sed -i'' "s:#!.*python:#!$PYTHON:" "$f"
  echo " by '$(head -n1 "$f")'"
done

########################################################################
# Write custom AppRun script
########################################################################

cp "$TRAVIS_BUILD_DIR/Documentation/static/logo.svg" "$APP.svg"

cat -- > "$APP.desktop" <<EOF
[Desktop Entry]
Version=1.0
Name=$APP
Icon=$APP
Comment=Medical Image Registration ToolKit ($VERSION)
Exec=$PYTHON usr/bin/mirtk
Terminal=true
Type=Application
Categories=Development
EOF

cat -- > "AppRun" <<EOF
#!/bin/bash
APPDIR="\$(dirname "\$(readlink -f "\$BASH_SOURCE")")"
export LD_LIBRARY_PATH="\$APPDIR/usr/lib:\$LD_LIBRARY_PATH"
$PYTHON "\$APPDIR/usr/bin/mirtk" "\$@"
EOF
chmod a+x AppRun

echo "Populating AppDir... done"
tree -a 2> /dev/null || true

########################################################################
# Now packaging AppDir as an AppImage
########################################################################

echo
echo "Generating AppImages..."
APPIMAGES=()

cd "$APPDIR/.."
mkdir -p "$OUTDIR"
wget "https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-${ARCH}.AppImage"
chmod a+x appimagetool-${ARCH}.AppImage

if [ $AppImage_VIEWER = ON ]; then
  echo "Generating AppImage including view binary..."
  APPIMAGE="$APP+view-$RELEASE-$ARCH"
  GLIBC_VERSION=$(glibc_needed)
  [ -z "$GLIBC_VERSION" ] || APPIMAGE="$APPIMAGE-glibc$GLIBC_VERSION"
  APPIMAGE="$APPIMAGE.AppImage"
  rm -f "$OUTDIR/$APPIMAGE" 2> /dev/null || true
  ./appimagetool-${ARCH}.AppImage -n -v "$APPDIR" "$OUTDIR/$APPIMAGE"
  APPIMAGES=("${APPIMAGES[@]}" "$APPIMAGE")
  echo "Generating AppImage including view binary... done"

  cd "$APPDIR/usr/lib"
  rm -f tools/view \
        libgsl* \
        libGL* \
        libX* \
        libxcb* \
        libxshmfence* \
        libfreetype* \
        libglapi* \
  rm -rf mesa || true
  cd "$APPDIR/.."
fi

echo "Generating AppImage excluding view binary..."
APPIMAGE="$APP-$RELEASE-$ARCH"
GLIBC_VERSION=$(glibc_needed)
[ -z "$GLIBC_VERSION" ] || APPIMAGE="$APPIMAGE-glibc$GLIBC_VERSION"
APPIMAGE="$APPIMAGE.AppImage"
rm -f "$OUTDIR/$APPIMAGE" 2> /dev/null || true
./appimagetool-${ARCH}.AppImage -n -v "$APPDIR" "$OUTDIR/$APPIMAGE"
APPIMAGES=("${APPIMAGES[@]}" "$APPIMAGE")
echo "Generating AppImage excluding view binary... done"

echo "Generating AppImages... done"

########################################################################
# Clean up
########################################################################

if [ $TRAVIS = OFF ]; then
  rm -rf "$APPDIR" "/tmp/$APP" || true
fi

########################################################################
# Upload AppImage
########################################################################

if [ $TRANSFER = ON ]; then
  echo
  echo "Transferring AppImages ..."
  for APPIMAGE in "${APPIMAGES[@]}"; do
    transfer "$OUTDIR/$APPIMAGE"
    echo
  done
  echo "Transferring AppImages... done"
fi

cat -- <<EOF

You can test the AppImage using a Docker container, e.g., a fresh Ubuntu installation.

# Without extracting the AppImage using FUSE
(Not recommended, see https://github.com/AppImage/AppImageKit/wiki/FUSE#docker)

Pull and run Docker container with necessary privileges for FUSE:
> docker run --rm -it --privileged=true --cap-add MKNOD --cap-add SYS_ADMIN --device /dev/fuse ubuntu

Inside the Docker container run:
> URL="Copy and paste MIRTK AppImage transfer.sh, bintray.com, or github.com URL here"
> apt-get update && apt-get install -y sshfs wget python
> wget -O mirtk "\$URL" && chmod a+x mirtk
> ./mirtk -h

# Extracting the AppImage and calling squashfs-root/AppRun

Pull and run Docker container:
> docker run --rm -it ubuntu

Inside the Docker container run:
> URL="Copy and paste MIRTK AppImage transfer.sh, bintray.com, or github.com URL here"
> apt-get update && apt-get install -y wget python
> wget -O MIRTK.AppImage "\$URL" && chmod a+x MIRTK.AppImage
> ./MIRTK.AppImage --appimage-extract
> ./squashfs-root/AppRun -h
EOF
