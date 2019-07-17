#!/usr/bin/env bash
set -vex

################
# DEPENDENCIES #
################

## Load modules
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module purge

module load gcc
module load ccache

module load meson
module load ninja

module load boost
module load hdf5-tools
module load cram

case "${bamboo_planRepository_branchName}" in
  master)
    module load libblasr/master
    module load pbbam/master
    ;;
  *)
    module load libblasr/develop
    module load pbbam/develop
    ;;
esac


BOOST_ROOT="${BOOST_ROOT%/include}"
# unset these variables to have meson discover all
# boost-dependent variables from BOOST_ROOT alone
unset BOOST_INCLUDEDIR
unset BOOST_LIBRARYDIR

export CC="ccache gcc"
export CXX="ccache g++"
export CCACHE_BASEDIR="${PWD}"

if [[ $USER == bamboo ]]; then
  export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.${bamboo_shortPlanKey}.ccachedir
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi

case "${bamboo_planRepository_branchName}" in
  develop|master)
    export PREFIX_ARG="/mnt/software/b/blasr/${bamboo_planRepository_branchName}"
    export BUILD_NUMBER="${bamboo_globalBuildNumber:-0}"
    ;;
  *)
    export BUILD_NUMBER="0"
    ;;
esac

# call the main build+test scripts
export CURRENT_BUILD_DIR="build"
export ENABLED_TESTS="true"

# TODO(dseifert)
# HDF5 doesn't have pkg-config files yet
export CPPFLAGS="${HDF5_CFLAGS}"
export LDFLAGS="-static-libstdc++ -static-libgcc ${HDF5_LIBS} -Wl,-rpath-link,${LD_LIBRARY_PATH}"

bash scripts/ci/build.sh
bash scripts/ci/test.sh

if [[ -z ${PREFIX_ARG+x} ]]; then
  echo "Not installing anything (branch: ${bamboo_planRepository_branchName}), exiting."
  exit 0
fi

bash scripts/ci/install.sh
