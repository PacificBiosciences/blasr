## Compile BLASR and util binaries using cmake

Make sure that you are using cmake >=3.7, have [ninja](https://ninja-build.org/)
installed, and always start from an empty build subdirectory!

    git clone git://github.com/PacificBiosciences/blasr.git && cd blasr
    git submodule update --init --remote
    mkdir build && cd build
    cmake -GNinja .. && ninja

Where is the binary?

    blasr/build/blasr

Can I install it system-wide?

    ninja install

Can I use make instead of ninja, even though it's officially not supported?
Yes, but make does not resolve dependencies correctly;
because of that, you need to run make twice.

    cmake .. && make -j && make -j

Can I install all binaries to a custom location?

    cmake -GNinja -DCMAKE_INSTALL_PREFIX=/my/path .. && ninja install

Is your HDF5 in a custom location?

    cmake -GNinja -DHDF5_ROOT=/your/location/hdf-1.8.16 ..

Are HDF$ libraries and include folders in different locations?

    cmake -GNinja -DHDF5_LIBRARIES=/your/location/hdf-1.8.16/lib
                  -DHDF5_INCLUDE_DIRS=/other/location/hdf-1.8.16/include ..

To change the location of htslib, add:
    -DHTSLIB_INCLUDE_DIR=/path/to/htslib -DHTSLIB_LIB_DIR=/path/to/htslib


Prefer a custom libz implementation?

    cmake -GNinja -DZLIB_INCLUDE_DIRS=/your/location/zlib/include \
                  -DZLIB_LIBRARIES=/your/location/zlib/libz.so ..

Are htslib libraries and include folders in different locations?
    cmake -GNinja -DHTSLIB_INCLUDE_DIR=/your/location/htslib \
                  -DHTSLIB_LIB_DIR=/your/location/htslib ..

Are you on MacOS? Install homebrew and install following dependencies and
use the default instructions above:

    brew install ninja hdf5 cmake
