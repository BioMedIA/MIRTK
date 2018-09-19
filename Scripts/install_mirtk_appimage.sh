#!/bin/bash

INSTALL_PREFIX="$HOME/opt/mirtk-appimage"
SUDO=''  # set to empty string if prefix is writable by $USER

APPIMAGE="$1"
if [ -z "$APPIMAGE" ]; then
  GLIBC_VERSION="$(ldd --version | head -n1 | rev | cut -d' ' -f1 | rev)"
  GLIBC_VERSION_MAJOR="${GLIBC_VERSION/.*}"
  GLIBC_VERSION_MINOR="${GLIBC_VERSION/*.}"
  if [ $GLIBC_VERSION_MAJOR -lt 2 ] || [ $GLIBC_VERSION_MAJOR -eq 2 -a $GLIBC_VERSION_MINOR -lt 14 ]; then
    cat -- 2>&1 <<EOF
Your glibc version is $GLIBC_VERSION, but at least version 2.14 is required!
Consider building MIRTK from source code or use the Docker container.
See also: http://mirtk.github.io/getstarted.html#install-the-software
EOF
    exit 1
  else
    APPIMAGE=MIRTK-latest-x86_64-glibc2.14.AppImage
  fi
else
  name="$(basename "$APPIMAGE")"
  if [ "${name:0:5}" != MIRTK -o "${name//*.}" != AppImage ]; then
    echo "Provided AppImage file must be named MIRTK*.AppImage!" 2>&1
    exit 1
  fi
  if [ ! -f "$APPIMAGE" -a "$APPIMAGE" != "$name" ]; then
    echo "Specified AppImage file does not exist!" 2>&1
    echo "To download an AppImage from bintray.com, specify file name only!" 2>&1
    exit 1
  fi
fi

if [ -f "$APPIMAGE" ]; then
  cp "$APPIMAGE" /tmp/mirtk
elif [ $(which wget 2> /dev/null) ]; then
  wget -O /tmp/mirtk https://bintray.com/schuhschuh/AppImages/download_file?file_path=$APPIMAGE
elif [ $(which curl 2> /dev/null) ]; then
  curl -o /tmp/mirtk https://bintray.com/schuhschuh/AppImages/download_file?file_path=$APPIMAGE
else
  cat -- <<EOF
Neither wget nor curl is installed!
Download the following file and call this script with downloaded file as argument:
https://bintray.com/schuhschuh/AppImages/download_file?file_path=MIRTK-latest-x86_64-glibc$GLIBC_MINIMUM.AppImage
EOF
fi

chmod a+x /tmp/mirtk && $SUDO mkdir -p "$INSTALL_PREFIX/bin" && $SUDO cp /tmp/mirtk "$INSTALL_PREFIX/bin/mirtk" && rm -f /tmp/mirtk
