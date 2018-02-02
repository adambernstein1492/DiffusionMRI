#!/bin/sh
set +x
cd ${1}/mcell-3
autoheader
aclocal
automake --add-missing
autoconf
echo $TARGET $CC $CFLAGS
if [ "x$TARGET" = "xmic" ]; then
  ./configure CC=icc CXX=icpc CFLAGS="-mkl -std=gnu99 -DDO_MPI -DMRI_DIFF_SIM -O3 -Wall" --prefix=${1}/mcell
  sed -i -e 's/-mkl/-mmic -mkl/g' Makefile
else
  if [ "x$CC" = "xicc" ]; then
    ./configure CC=icc CXX=icpc CFLAGS="-mkl -std=gnu99 -DDO_MPI -DMRI_DIFF_SIM -O3 -Wall" --prefix=${1}/mcell
  else
    ./configure CC=gcc CXX=g++ CFLAGS="-std=gnu99 -DDO_MPI -DMRI_DIFF_SIM -O3 -Wall" --prefix=${1}/mcell
  fi
fi
